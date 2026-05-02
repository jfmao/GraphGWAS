"""Parallel computing infrastructure for GraphGWAS.

Provides process-level parallelism for chromosome-level fan-out (genome-wide scans),
window-level fan-out (epistasis), and permutation parallelism.

Core challenge: neo4j.Driver is thread-safe but not picklable across process boundaries.
Solution: ConnectionParams dataclass carries credentials; each worker creates its own
connection inside the subprocess.

Usage:
    from graphgwas.parallel import parallel_genome_scan, connection_params_from_conn
    results = parallel_genome_scan(conn, chromosomes, method="auto", max_workers=8)
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# Connection parameter transport (picklable across process boundaries)
# ---------------------------------------------------------------------------

@dataclass
class ConnectionParams:
    """Picklable Neo4j connection parameters for worker processes."""
    uri: str
    user: str
    password: str
    database: str | None


def connection_params_from_conn(conn) -> ConnectionParams:
    """Extract picklable ConnectionParams from a GraphGWASConnection."""
    return ConnectionParams(
        uri=conn._uri,
        user=conn._user,
        password=conn._password,
        database=conn._database,
    )


def _default_workers(task_count: int, cap: int = 8) -> int:
    """Determine default worker count: min(task_count, cpu_count, cap)."""
    cpu = os.cpu_count() or 4
    return max(1, min(task_count, cpu, cap))


# ---------------------------------------------------------------------------
# Worker functions (module-level for picklability)
# ---------------------------------------------------------------------------

def _scan_chromosome_worker(params: ConnectionParams, chr: str,
                            method: str,
                            covariate_names: list[str] | None,
                            worker_index: int = 0) -> list[dict]:
    """Worker: run single_locus_scan on one chromosome."""
    import time
    from .db import GraphGWASConnection
    from .assoc import single_locus_scan

    # Stagger connection creation to avoid thundering herd
    if worker_index > 0:
        time.sleep(worker_index * 0.2)

    with GraphGWASConnection(params.uri, params.user,
                             params.password, params.database) as conn:
        return single_locus_scan(conn, chr, method=method,
                                 covariate_names=covariate_names, verbose=False)


def _epistasis_window_worker(params: ConnectionParams, chr: str,
                             w_start: int, w_end: int,
                             min_cocarriers: int, min_enrichment: float,
                             n_permutations: int) -> list[dict]:
    """Worker: process one epistasis window (build graph + detect + test)."""
    from .db import GraphGWASConnection
    from .epistasis import (
        build_cooccurrence_graph,
        detect_epistatic_modules,
        test_epistatic_module,
    )

    with GraphGWASConnection(params.uri, params.user,
                             params.password, params.database) as conn:
        G = build_cooccurrence_graph(conn, chr, w_start, w_end,
                                     min_cocarriers=min_cocarriers,
                                     min_enrichment=min_enrichment,
                                     verbose=False)
        if G.number_of_nodes() < 2:
            return []

        modules = detect_epistatic_modules(G)
        results = []
        for module in modules:
            result = test_epistatic_module(conn, module,
                                          n_permutations=n_permutations,
                                          verbose=False)
            result["chr"] = chr
            result["window_start"] = w_start
            result["window_end"] = w_end
            results.append(result)
        return results


def _permutation_batch_worker(all_carriers_list: list[np.ndarray],
                              n_case: int, n_ctrl: int,
                              seed: int, n_perms: int) -> np.ndarray:
    """Worker: run a batch of permutations for epistasis module testing.

    Pure compute — no database access needed. Carrier arrays passed directly.
    Returns array of null module scores.
    """
    from .epistasis import compute_module_score

    rng = np.random.default_rng(seed)
    n_all = n_case + n_ctrl
    null_scores = np.empty(n_perms)

    for i in range(n_perms):
        perm_labels = rng.permutation(n_all)
        perm_case_mask = perm_labels < n_case

        perm_case_carriers = [cs[perm_case_mask] for cs in all_carriers_list]
        perm_ctrl_carriers = [cs[~perm_case_mask] for cs in all_carriers_list]

        null_scores[i] = compute_module_score(
            perm_case_carriers, perm_ctrl_carriers, n_case, n_ctrl
        )

    return null_scores


def _flow_permutation_worker(all_sample_nodes: list[str],
                             variant_gene_edges: list[dict],
                             n_case: int, seed: int,
                             n_perms: int) -> dict[str, list[float]]:
    """Worker: run a batch of flow permutations.

    Rebuilds Source edges with permuted case labels and computes pathway flows.
    Takes serializable network data, not nx.DiGraph (not picklable reliably).
    """
    from .flow import compute_pathway_flows

    rng = np.random.default_rng(seed)
    null_flows: dict[str, list[float]] = {}

    for _ in range(n_perms):
        # Rebuild network from serialized data with permuted labels
        G = _rebuild_flow_network(variant_gene_edges, all_sample_nodes,
                                  set(rng.choice(all_sample_nodes, size=n_case,
                                                 replace=False)))

        perm_results = compute_pathway_flows(G, verbose=False)
        for r in perm_results:
            pw = r["pathway"]
            if pw not in null_flows:
                null_flows[pw] = []
            null_flows[pw].append(r["flow_value"])

    return null_flows


def _rebuild_flow_network(edge_data: list[dict], all_sample_nodes: list[str],
                          case_nodes: set[str]) -> "nx.DiGraph":
    """Rebuild a flow network from serialized edge data with given case labels."""
    import networkx as nx

    G = nx.DiGraph()
    G.add_node("SOURCE", type="source")

    for e in edge_data:
        G.add_node(e["src"], type=e["src_type"])
        G.add_node(e["dst"], type=e["dst_type"])
        G.add_edge(e["src"], e["dst"], capacity=e["capacity"])

    # Source → case sample edges
    for sn in all_sample_nodes:
        if sn not in G:
            G.add_node(sn, type="sample")
        if sn in case_nodes:
            G.add_edge("SOURCE", sn, capacity=1.0)

    return G


# ---------------------------------------------------------------------------
# Orchestrators
# ---------------------------------------------------------------------------

def parallel_genome_scan(conn, chromosomes: list[str],
                         method: str = "auto",
                         covariate_names: list[str] | None = None,
                         max_workers: int | None = None,
                         verbose: bool = True) -> list[dict]:
    """Fan out single_locus_scan across chromosomes using ProcessPoolExecutor.

    Each worker creates its own Neo4j connection.
    Results are collected and merged.
    """
    params = connection_params_from_conn(conn)
    n_workers = max_workers or _default_workers(len(chromosomes))

    if verbose:
        print(f"Parallel genome scan: {len(chromosomes)} chromosomes, "
              f"{n_workers} workers")

    all_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(_scan_chromosome_worker, params, chr, method,
                        covariate_names, i): chr
            for i, chr in enumerate(chromosomes)
        }
        for future in as_completed(futures):
            chr_name = futures[future]
            try:
                results = future.result()
                all_results.extend(results)
                if verbose:
                    sig = sum(1 for r in results if r["p_value"] < 5e-8)
                    print(f"  {chr_name}: {len(results)} variants, "
                          f"{sig} significant")
            except Exception as e:
                if verbose:
                    print(f"  {chr_name}: FAILED — {e}")

    return all_results


def parallel_epistasis_scan(conn, chr: str,
                            windows: list[tuple[int, int]],
                            min_cocarriers: int = 3,
                            min_enrichment: float = 1.5,
                            n_permutations: int = 1000,
                            max_workers: int | None = None,
                            verbose: bool = True) -> list[dict]:
    """Fan out epistasis scan across genomic windows."""
    params = connection_params_from_conn(conn)
    n_workers = max_workers or _default_workers(len(windows))

    if verbose:
        print(f"Parallel epistasis scan: {len(windows)} windows, "
              f"{n_workers} workers")

    all_modules = []
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(_epistasis_window_worker, params, chr,
                        ws, we, min_cocarriers, min_enrichment,
                        n_permutations): (ws, we)
            for ws, we in windows
        }
        for future in as_completed(futures):
            window = futures[future]
            try:
                modules = future.result()
                all_modules.extend(modules)
                if verbose and modules:
                    print(f"  Window {window[0]}-{window[1]}: "
                          f"{len(modules)} modules")
            except Exception as e:
                if verbose:
                    print(f"  Window {window[0]}-{window[1]}: FAILED — {e}")

    all_modules.sort(key=lambda m: m["p_value"])
    return all_modules


def parallel_permutation_test(all_carriers_list: list[np.ndarray],
                              n_case: int, n_ctrl: int,
                              observed_score: float,
                              n_permutations: int = 1000,
                              max_workers: int | None = None) -> float:
    """Parallelize permutation test for epistasis module significance.

    Pure compute — no database access. Splits permutations across workers.
    Returns empirical p-value.
    """
    n_workers = max_workers or _default_workers(n_permutations, cap=8)
    perms_per_worker = n_permutations // n_workers
    remainder = n_permutations % n_workers

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = []
        for i in range(n_workers):
            n = perms_per_worker + (1 if i < remainder else 0)
            if n > 0:
                futures.append(
                    pool.submit(_permutation_batch_worker,
                                all_carriers_list, n_case, n_ctrl,
                                seed=42 + i, n_perms=n)
                )

        null_scores = np.concatenate([f.result() for f in futures])

    p_value = float((np.sum(null_scores >= observed_score) + 1) /
                    (len(null_scores) + 1))
    return p_value


def serialize_flow_network(G) -> tuple[list[dict], list[str]]:
    """Serialize a NetworkX DiGraph's edges for cross-process transport.

    Returns (edge_data, sample_node_list).
    """
    edge_data = []
    sample_nodes = []

    for u, v, data in G.edges(data=True):
        if u == "SOURCE":
            continue  # Source edges are rebuilt per permutation
        src_type = G.nodes[u].get("type", "unknown")
        dst_type = G.nodes[v].get("type", "unknown")
        edge_data.append({
            "src": u, "dst": v,
            "src_type": src_type, "dst_type": dst_type,
            "capacity": data.get("capacity", 1.0),
        })

    for n, d in G.nodes(data=True):
        if d.get("type") == "sample":
            sample_nodes.append(n)

    return edge_data, sample_nodes


def parallel_flow_permutations(flow_network,
                               n_case: int,
                               observed_flows: dict[str, float],
                               n_permutations: int = 100,
                               max_workers: int | None = None,
                               verbose: bool = True) -> dict[str, float]:
    """Parallelize flow permutation test across workers.

    Returns dict of pathway → p-value.
    """
    edge_data, sample_nodes = serialize_flow_network(flow_network)
    n_workers = max_workers or _default_workers(n_permutations, cap=4)
    perms_per_worker = n_permutations // n_workers
    remainder = n_permutations % n_workers

    if verbose:
        print(f"Parallel flow permutations: {n_permutations} perms, "
              f"{n_workers} workers")

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = []
        for i in range(n_workers):
            n = perms_per_worker + (1 if i < remainder else 0)
            if n > 0:
                futures.append(
                    pool.submit(_flow_permutation_worker,
                                sample_nodes, edge_data, n_case,
                                seed=42 + i, n_perms=n)
                )

        # Merge null distributions from all workers
        merged_nulls: dict[str, list[float]] = {}
        for f in futures:
            worker_nulls = f.result()
            for pw, vals in worker_nulls.items():
                if pw not in merged_nulls:
                    merged_nulls[pw] = []
                merged_nulls[pw].extend(vals)

    # Compute p-values
    p_values = {}
    for pw, obs_flow in observed_flows.items():
        null_arr = np.array(merged_nulls.get(pw, [0.0]))
        p_values[pw] = float(
            (np.sum(null_arr >= obs_flow) + 1) / (len(null_arr) + 1)
        )

    return p_values
