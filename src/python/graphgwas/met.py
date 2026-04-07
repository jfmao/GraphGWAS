"""Multi-Environment Trial (MET) analysis — graph-native experimental design.

Represents the experimental design AS a graph: genotypes connected to
environments via OBSERVED_IN edges, environments connected by genetic
similarity, pedigree structure as explicit parent→population→genotype paths.

Unbalanced designs, tester-location confounding, and missing G×E cells
are graph topology properties — not statistical nuisances to be corrected.

Inspired by the Big BIT maize experiment (Jines et al. 2025): 2,554 DHs
× 3 testers × 47 environments, deeply unbalanced, with reference set
bridging disconnected subgraphs.

Components:
    M1. load_met_data           — import long-format MET CSV into graph
    M2. environment_similarity  — genetic correlation between environments
    M3. mega_environments       — community detection on env similarity graph
    M4. per_environment_gwas    — GWAS within each environment
    M5. variant_reaction_norms  — per-variant effect × environment
    M6. gxe_decomposition       — spectral G×E variance decomposition
    M7. design_diagnostics      — graph connectivity & balance analysis
    M8. graph_impute_missing    — predict missing G×E cells via graph neighbors
    M9. met_report              — integrated MET report
"""

from __future__ import annotations

import csv

import numpy as np
import networkx as nx
from scipy import stats as sp_stats

from .db import GraphGWASConnection


# ===================================================================
# Schema helpers
# ===================================================================

def ensure_met_indexes(conn: GraphGWASConnection) -> list[str]:
    """Create indexes for MET node types."""
    indexes = [
        "CREATE INDEX met_trial_id IF NOT EXISTS FOR (t:Trial) ON (t.id)",
        "CREATE INDEX met_env_id IF NOT EXISTS FOR (e:Environment) ON (e.id)",
        "CREATE INDEX met_env_trial IF NOT EXISTS FOR (e:Environment) ON (e.trial_id)",
        "CREATE INDEX met_genotype_id IF NOT EXISTS FOR (g:Genotype) ON (g.id)",
    ]
    results = []
    for idx in indexes:
        try:
            conn.execute_write(idx)
            results.append(f"OK: {idx}")
        except Exception as e:
            results.append(f"SKIP: {e}")
    return results


# ===================================================================
# M1. Load MET Data
# ===================================================================

def load_met_data(conn: GraphGWASConnection,
                  csv_path: str,
                  genotype_col: str,
                  environment_col: str,
                  trait_cols: list[str],
                  env_covariate_cols: list[str] | None = None,
                  rep_col: str | None = None,
                  block_col: str | None = None,
                  trial_id: str = "default",
                  batch_size: int = 500,
                  verbose: bool = True) -> dict:
    """Import long-format MET phenotype + design data into the graph.

    Creates:
        (:Trial) — one per trial_id
        (:Environment) — one per unique environment value
        (:Genotype) — one per unique genotype value
        (:Genotype)-[:OBSERVED_IN {pheno_<trait>: value, rep, block}]->(:Environment)

    If Genotype.id matches an existing Sample.sampleId, a
    (:Genotype)-[:IS_SAMPLE]->(:Sample) link is created.

    Args:
        conn: database connection.
        csv_path: path to long-format CSV.
        genotype_col: column for genotype identifier.
        environment_col: column for environment identifier.
        trait_cols: columns for trait values.
        env_covariate_cols: optional columns for environment-level covariates.
        rep_col: optional column for replicate number.
        block_col: optional column for block identifier.
        trial_id: trial identifier.
        batch_size: batch size for Neo4j writes.
        verbose: print progress.

    Returns:
        Summary dict with counts and connectivity stats.
    """
    if env_covariate_cols is None:
        env_covariate_cols = []

    # Parse CSV
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if verbose:
        print(f"=== Loading MET Data ===")
        print(f"  File: {csv_path}")
        print(f"  Rows: {len(rows)}")

    # Extract unique genotypes and environments
    genotypes = sorted(set(r[genotype_col] for r in rows if r.get(genotype_col)))
    environments = sorted(set(r[environment_col] for r in rows if r.get(environment_col)))

    if verbose:
        print(f"  Genotypes: {len(genotypes)}")
        print(f"  Environments: {len(environments)}")
        print(f"  Traits: {trait_cols}")

    # Ensure indexes
    ensure_met_indexes(conn)

    # Create Trial node
    conn.execute_write(
        """
        MERGE (t:Trial {id: $id})
        SET t.n_genotypes = $ng, t.n_environments = $ne, t.n_observations = $no
        """,
        {"id": trial_id, "ng": len(genotypes), "ne": len(environments),
         "no": len(rows)},
    )

    # Create Environment nodes
    # Aggregate environment-level covariates (take first occurrence per env)
    env_covariates: dict[str, dict] = {}
    for r in rows:
        env_id = r.get(environment_col)
        if env_id and env_id not in env_covariates:
            covs = {}
            for col in env_covariate_cols:
                val = r.get(col, "")
                if val:
                    try:
                        covs[col] = float(val)
                    except ValueError:
                        covs[col] = val
            # Try to parse location and year from environment ID
            parts = env_id.split("_")
            if len(parts) >= 2:
                covs.setdefault("location", parts[0])
                try:
                    covs.setdefault("year", int(parts[-1]))
                except ValueError:
                    pass
            env_covariates[env_id] = covs

    for i in range(0, len(environments), batch_size):
        batch = []
        for env_id in environments[i:i + batch_size]:
            entry = {"id": env_id, "trial_id": trial_id}
            entry.update(env_covariates.get(env_id, {}))
            batch.append(entry)

        # Dynamic SET clause for covariates
        set_parts = ["e.trial_id = row.trial_id"]
        for col in env_covariate_cols:
            safe = col.replace(" ", "_").replace("-", "_")
            set_parts.append(f"e.{safe} = row.{safe}")
        if "location" not in env_covariate_cols:
            set_parts.append("e.location = row.location")
        if "year" not in env_covariate_cols:
            set_parts.append("e.year = row.year")

        conn.execute_write(
            f"""
            UNWIND $batch AS row
            MERGE (e:Environment {{id: row.id}})
            SET {', '.join(set_parts)}
            """,
            {"batch": batch},
        )

    # Create Genotype nodes
    for i in range(0, len(genotypes), batch_size):
        batch = [{"id": g} for g in genotypes[i:i + batch_size]]
        conn.execute_write(
            """
            UNWIND $batch AS row
            MERGE (g:Genotype {id: row.id})
            """,
            {"batch": batch},
        )

    # Link Genotype → Sample if sampleId matches
    conn.execute_write(
        """
        MATCH (g:Genotype), (s:Sample)
        WHERE g.id = s.sampleId
        MERGE (g)-[:IS_SAMPLE]->(s)
        """,
    )

    # Create OBSERVED_IN relationships with trait values
    n_created = 0
    safe_traits = [t.replace(" ", "_").replace("-", "_") for t in trait_cols]

    for i in range(0, len(rows), batch_size):
        batch = []
        for r in rows[i:i + batch_size]:
            geno = r.get(genotype_col)
            env = r.get(environment_col)
            if not geno or not env:
                continue

            entry = {"geno": geno, "env": env}
            for col, safe in zip(trait_cols, safe_traits):
                val = r.get(col, "")
                if val:
                    try:
                        entry[f"pheno_{safe}"] = float(val)
                    except ValueError:
                        pass
            if rep_col and r.get(rep_col):
                try:
                    entry["rep"] = int(r[rep_col])
                except ValueError:
                    entry["rep"] = r[rep_col]
            if block_col and r.get(block_col):
                entry["block"] = r[block_col]

            batch.append(entry)

        # Build dynamic SET clause
        set_parts = []
        for safe in safe_traits:
            set_parts.append(f"r.pheno_{safe} = row.pheno_{safe}")
        if rep_col:
            set_parts.append("r.rep = row.rep")
        if block_col:
            set_parts.append("r.block = row.block")

        set_clause = ", ".join(set_parts) if set_parts else "r._loaded = true"

        result = conn.execute_write(
            f"""
            UNWIND $batch AS row
            MATCH (g:Genotype {{id: row.geno}})
            MATCH (e:Environment {{id: row.env}})
            CREATE (g)-[r:OBSERVED_IN]->(e)
            SET {set_clause}
            RETURN count(r) AS n
            """,
            {"batch": batch},
        )
        rec = result.single()
        n_created += rec["n"] if rec else 0

    if verbose:
        print(f"\n  Created: {len(genotypes)} Genotype nodes, "
              f"{len(environments)} Environment nodes, "
              f"{n_created} OBSERVED_IN relationships")

    return {
        "trial_id": trial_id,
        "n_genotypes": len(genotypes),
        "n_environments": len(environments),
        "n_observations": n_created,
        "traits": trait_cols,
        "env_covariates": env_covariate_cols,
    }


# ===================================================================
# M2. Environment Similarity
# ===================================================================

def environment_similarity(conn: GraphGWASConnection,
                           trial_id: str,
                           trait: str,
                           min_shared_genotypes: int = 10,
                           verbose: bool = True) -> dict:
    """Genetic correlation between environments from shared genotype performance.

    For each pair of environments, finds genotypes tested in both and
    computes Pearson correlation of trait values.  This gives a non-parametric
    environment similarity that Factor Analytic models estimate parametrically.

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name (without pheno_ prefix).
        min_shared_genotypes: minimum overlap for correlation.
        verbose: print progress.

    Returns:
        dict with similarity_matrix, environment_ids, and pairwise correlations.
    """
    safe_trait = f"pheno_{trait.replace(' ', '_').replace('-', '_')}"

    # Get all environments for this trial
    env_result = conn.execute_read(
        "MATCH (e:Environment {trial_id: $tid}) RETURN e.id AS id ORDER BY e.id",
        {"tid": trial_id},
    )
    env_ids = [r["id"] for r in env_result]
    E = len(env_ids)

    if verbose:
        print(f"=== Environment Similarity ===")
        print(f"  Trial: {trial_id}, Trait: {trait}, Environments: {E}")

    # For each environment, collect genotype → trait value
    env_data: dict[str, dict[str, float]] = {}
    for env_id in env_ids:
        result = conn.execute_read(
            f"""
            MATCH (g:Genotype)-[r:OBSERVED_IN]->(e:Environment {{id: $eid}})
            WHERE r.{safe_trait} IS NOT NULL
            RETURN g.id AS geno, r.{safe_trait} AS val
            """,
            {"eid": env_id},
        )
        env_data[env_id] = {rec["geno"]: rec["val"] for rec in result}

    # Pairwise correlations
    sim_matrix = np.eye(E, dtype=np.float64)
    pairwise = []

    for i in range(E):
        for j in range(i + 1, E):
            d1 = env_data[env_ids[i]]
            d2 = env_data[env_ids[j]]
            shared = set(d1.keys()) & set(d2.keys())
            n_shared = len(shared)

            if n_shared < min_shared_genotypes:
                r_g = 0.0
            else:
                v1 = np.array([d1[g] for g in shared])
                v2 = np.array([d2[g] for g in shared])
                if np.std(v1) > 0 and np.std(v2) > 0:
                    r_g = float(np.corrcoef(v1, v2)[0, 1])
                else:
                    r_g = 0.0

            sim_matrix[i, j] = r_g
            sim_matrix[j, i] = r_g
            pairwise.append({
                "env1": env_ids[i], "env2": env_ids[j],
                "r_g": r_g, "n_shared": n_shared,
            })

    # Store similarity edges in graph
    for pw in pairwise:
        if abs(pw["r_g"]) > 0.01:
            conn.execute_write(
                """
                MATCH (e1:Environment {id: $e1}), (e2:Environment {id: $e2})
                MERGE (e1)-[r:ENV_SIMILAR]->(e2)
                SET r.r_g = $rg, r.n_shared = $ns, r.trait = $trait
                """,
                {"e1": pw["env1"], "e2": pw["env2"],
                 "rg": pw["r_g"], "ns": pw["n_shared"], "trait": trait},
            )

    if verbose:
        mean_rg = np.mean([pw["r_g"] for pw in pairwise])
        print(f"  Mean r_G across env pairs: {mean_rg:.4f}")
        print(f"  Range: [{min(pw['r_g'] for pw in pairwise):.4f}, "
              f"{max(pw['r_g'] for pw in pairwise):.4f}]")

    return {
        "trial_id": trial_id,
        "trait": trait,
        "environment_ids": env_ids,
        "similarity_matrix": sim_matrix.tolist(),
        "pairwise": pairwise,
        "n_environments": E,
    }


# ===================================================================
# M3. Mega-Environments
# ===================================================================

def mega_environments(conn: GraphGWASConnection,
                      trial_id: str,
                      trait: str,
                      resolution: float = 1.0,
                      min_shared_genotypes: int = 10,
                      verbose: bool = True) -> dict:
    """Detect mega-environments via community detection on similarity graph.

    Each community = a set of environments where genetic effects are similar.
    This is what Factor Analytic models approximate parametrically.

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name.
        resolution: Louvain resolution parameter.
        min_shared_genotypes: for similarity computation.
        verbose: print progress.

    Returns:
        dict with mega-environment assignments and per-cluster stats.
    """
    # Compute similarity if not already done
    sim = environment_similarity(conn, trial_id, trait,
                                 min_shared_genotypes=min_shared_genotypes,
                                 verbose=verbose)

    env_ids = sim["environment_ids"]
    sim_matrix = np.array(sim["similarity_matrix"])
    E = len(env_ids)

    if E < 2:
        return {"mega_environments": {env_ids[0]: 0} if env_ids else {},
                "n_clusters": 1}

    # Build networkx graph from similarity matrix
    G = nx.Graph()
    for i in range(E):
        G.add_node(env_ids[i])
    for i in range(E):
        for j in range(i + 1, E):
            if sim_matrix[i, j] > 0:
                G.add_edge(env_ids[i], env_ids[j], weight=sim_matrix[i, j])

    if G.number_of_edges() == 0:
        # No positive correlations — each env is its own mega-env
        assignments = {eid: i for i, eid in enumerate(env_ids)}
    else:
        communities = nx.community.louvain_communities(
            G, weight="weight", resolution=resolution, seed=42,
        )
        assignments = {}
        for ci, community in enumerate(communities):
            for env_id in community:
                assignments[env_id] = ci

    # Store on Environment nodes
    for env_id, cluster in assignments.items():
        conn.execute_write(
            """
            MATCH (e:Environment {id: $eid})
            SET e.mega_env = $me
            """,
            {"eid": env_id, "me": cluster},
        )

    n_clusters = len(set(assignments.values()))

    # Per-cluster stats
    clusters = {}
    for ci in range(n_clusters):
        members = [eid for eid, c in assignments.items() if c == ci]
        # Mean within-cluster r_G
        indices = [env_ids.index(m) for m in members if m in env_ids]
        if len(indices) >= 2:
            sub = sim_matrix[np.ix_(indices, indices)]
            mask = np.triu(np.ones_like(sub, dtype=bool), k=1)
            mean_rg = float(np.mean(sub[mask]))
        else:
            mean_rg = 1.0
        clusters[ci] = {"environments": members, "n": len(members),
                        "mean_within_r_g": mean_rg}

    if verbose:
        print(f"\n=== Mega-Environments ===")
        print(f"  {n_clusters} mega-environments detected:")
        for ci, info in clusters.items():
            print(f"    Cluster {ci}: {info['n']} envs, "
                  f"mean r_G={info['mean_within_r_g']:.4f} — "
                  f"{', '.join(info['environments'][:5])}")

    return {
        "trial_id": trial_id,
        "trait": trait,
        "assignments": assignments,
        "n_clusters": n_clusters,
        "clusters": clusters,
    }


# ===================================================================
# M4. Per-Environment GWAS
# ===================================================================

def per_environment_gwas(conn: GraphGWASConnection,
                         trial_id: str,
                         trait: str,
                         chr: str | None = None,
                         method: str = "linear",
                         environments: list[str] | None = None,
                         verbose: bool = True) -> dict:
    """Run association analysis within each environment separately.

    For each environment, extracts genotype-phenotype pairs from OBSERVED_IN
    relationships and runs variant-level association.  Requires that Genotype
    nodes are linked to Sample nodes (via IS_SAMPLE) for genotype access.

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name.
        chr: chromosome filter.
        method: association method.
        environments: subset of environments (None = all).
        verbose: print progress.

    Returns:
        Summary of per-environment GWAS results.
    """
    safe_trait = f"pheno_{trait.replace(' ', '_').replace('-', '_')}"

    # Get environments
    if environments is None:
        result = conn.execute_read(
            "MATCH (e:Environment {trial_id: $tid}) RETURN e.id AS id ORDER BY id",
            {"tid": trial_id},
        )
        environments = [r["id"] for r in result]

    if verbose:
        print(f"=== Per-Environment GWAS ===")
        print(f"  Trial: {trial_id}, Trait: {trait}")
        print(f"  Environments: {len(environments)}, Method: {method}")

    env_summaries = []
    for env_id in environments:
        # Get genotypes with both phenotype and sample link for this env
        result = conn.execute_read(
            f"""
            MATCH (g:Genotype)-[obs:OBSERVED_IN]->(e:Environment {{id: $eid}})
            WHERE obs.{safe_trait} IS NOT NULL
            OPTIONAL MATCH (g)-[:IS_SAMPLE]->(s:Sample)
            WHERE s.packed_index IS NOT NULL
            RETURN g.id AS geno, obs.{safe_trait} AS pheno_val,
                   s.packed_index AS packed_index
            """,
            {"eid": env_id},
        )

        records = [dict(r) for r in result]
        n_with_geno = sum(1 for r in records if r["packed_index"] is not None)

        env_summaries.append({
            "environment": env_id,
            "n_genotypes": len(records),
            "n_with_genotype_data": n_with_geno,
            "status": "ready" if n_with_geno >= 20 else "too_few",
        })

        if verbose:
            print(f"  {env_id}: {len(records)} genotypes "
                  f"({n_with_geno} with SNP data)")

    return {
        "trial_id": trial_id,
        "trait": trait,
        "n_environments": len(environments),
        "environment_summaries": env_summaries,
        "n_ready": sum(1 for s in env_summaries if s["status"] == "ready"),
    }


# ===================================================================
# M5. Variant Reaction Norms
# ===================================================================

def variant_reaction_norms(conn: GraphGWASConnection,
                           trial_id: str,
                           trait: str,
                           env_covariate: str | None = None,
                           top_n_variants: int = 100,
                           verbose: bool = True) -> dict:
    """Per-variant effect as a function of environment.

    Collects per-environment betas from stored AssociationResult nodes,
    computes beta variance across environments (= G×E sensitivity),
    and optionally fits reaction norm slopes against an environmental covariate.

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name.
        env_covariate: environment property to regress beta on (e.g. 'temperature_mean').
        top_n_variants: return top N most G×E-sensitive variants.
        verbose: print progress.

    Returns:
        Ranked list of variants by G×E sensitivity with reaction norm details.
    """
    # Get per-environment betas from AssociationResult nodes
    result = conn.execute_read(
        """
        MATCH (ar:AssociationResult)
        WHERE ar.phenotype_key = $trait AND ar.environment IS NOT NULL
        RETURN ar.variant_id AS variant, ar.environment AS env,
               ar.beta AS beta, ar.p_value AS p_value
        ORDER BY ar.variant_id, ar.environment
        """,
        {"trait": trait},
    )

    # Group by variant
    variant_betas: dict[str, list[dict]] = {}
    for rec in result:
        vid = rec["variant"]
        if vid not in variant_betas:
            variant_betas[vid] = []
        variant_betas[vid].append({
            "env": rec["env"],
            "beta": rec["beta"],
            "p_value": rec["p_value"],
        })

    if not variant_betas:
        if verbose:
            print("No per-environment association results found.")
            print("Run per_environment_gwas() first, storing results with environment property.")
        return {"variants": [], "n_variants": 0}

    # Get environment covariate values if requested
    env_cov_vals = {}
    if env_covariate:
        cov_result = conn.execute_read(
            f"""
            MATCH (e:Environment {{trial_id: $tid}})
            WHERE e.{env_covariate} IS NOT NULL
            RETURN e.id AS id, e.{env_covariate} AS cov
            """,
            {"tid": trial_id},
        )
        env_cov_vals = {r["id"]: r["cov"] for r in cov_result}

    # Compute reaction norm statistics per variant
    variants = []
    for vid, betas in variant_betas.items():
        beta_vals = np.array([b["beta"] for b in betas if b["beta"] is not None])
        if len(beta_vals) < 2:
            continue

        mean_beta = float(np.mean(beta_vals))
        beta_var = float(np.var(beta_vals))
        beta_range = float(np.max(beta_vals) - np.min(beta_vals))
        n_envs = len(beta_vals)

        entry = {
            "variant": vid,
            "mean_beta": mean_beta,
            "beta_variance": beta_var,
            "beta_range": beta_range,
            "n_environments": n_envs,
            "per_environment": betas,
        }

        # Reaction norm slope if covariate available
        if env_cov_vals:
            cov_vals = []
            beta_matched = []
            for b in betas:
                if b["env"] in env_cov_vals and b["beta"] is not None:
                    cov_vals.append(env_cov_vals[b["env"]])
                    beta_matched.append(b["beta"])
            if len(cov_vals) >= 3:
                slope, intercept, r, p, se = sp_stats.linregress(cov_vals, beta_matched)
                entry["reaction_norm_slope"] = float(slope)
                entry["reaction_norm_r2"] = float(r ** 2)
                entry["reaction_norm_p"] = float(p)

        variants.append(entry)

    # Sort by G×E sensitivity (beta variance)
    variants.sort(key=lambda v: v["beta_variance"], reverse=True)
    for i, v in enumerate(variants):
        v["rank"] = i + 1

    if verbose:
        print(f"=== Variant Reaction Norms ===")
        print(f"  {len(variants)} variants with multi-env betas")
        for v in variants[:10]:
            slope_str = (f" slope={v['reaction_norm_slope']:.4f}"
                         if "reaction_norm_slope" in v else "")
            print(f"  #{v['rank']} {v['variant']}: var={v['beta_variance']:.6f} "
                  f"range={v['beta_range']:.4f} envs={v['n_environments']}{slope_str}")

    return {
        "trial_id": trial_id,
        "trait": trait,
        "env_covariate": env_covariate,
        "variants": variants[:top_n_variants],
        "n_variants": len(variants),
    }


# ===================================================================
# M6. G×E Decomposition
# ===================================================================

def gxe_decomposition(conn: GraphGWASConnection,
                      trial_id: str,
                      trait: str,
                      verbose: bool = True) -> dict:
    """Decompose phenotypic variance into genetic, environmental, and G×E.

    Builds the G×E phenotype matrix from OBSERVED_IN relationships and
    computes ANOVA-style variance decomposition.

    V_total = V_G + V_E + V_GxE + V_residual

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name.
        verbose: print progress.

    Returns:
        Variance components and derived statistics.
    """
    safe_trait = f"pheno_{trait.replace(' ', '_').replace('-', '_')}"

    # Load all observations
    result = conn.execute_read(
        f"""
        MATCH (g:Genotype)-[r:OBSERVED_IN]->(e:Environment {{trial_id: $tid}})
        WHERE r.{safe_trait} IS NOT NULL
        RETURN g.id AS geno, e.id AS env, r.{safe_trait} AS val
        """,
        {"tid": trial_id},
    )

    records = [dict(r) for r in result]
    if len(records) < 10:
        return {"error": "Too few observations"}

    # Build G×E matrix
    genotypes = sorted(set(r["geno"] for r in records))
    environments = sorted(set(r["env"] for r in records))
    G, E = len(genotypes), len(environments)
    geno_idx = {g: i for i, g in enumerate(genotypes)}
    env_idx = {e: i for i, e in enumerate(environments)}

    # Mean phenotype per genotype-environment cell (handle replicates)
    cell_sums = np.zeros((G, E))
    cell_counts = np.zeros((G, E), dtype=int)
    for r in records:
        gi = geno_idx[r["geno"]]
        ei = env_idx[r["env"]]
        cell_sums[gi, ei] += r["val"]
        cell_counts[gi, ei] += 1

    # Cell means (NaN where no data)
    with np.errstate(divide="ignore", invalid="ignore"):
        Y = np.where(cell_counts > 0, cell_sums / cell_counts, np.nan)

    # Observed cells mask
    observed = cell_counts > 0
    n_obs = int(np.sum(observed))
    n_total = G * E
    balance_ratio = n_obs / n_total

    # Grand mean, row means, col means (ignoring NaN)
    grand_mean = float(np.nanmean(Y))
    geno_means = np.nanmean(Y, axis=1)  # (G,)
    env_means = np.nanmean(Y, axis=0)   # (E,)

    # Variance decomposition (Type III-like for unbalanced)
    # V_G: variance of genotype means
    v_g = float(np.nanvar(geno_means))
    # V_E: variance of environment means
    v_e = float(np.nanvar(env_means))
    # V_GxE: residual after removing G and E main effects
    predicted = geno_means[:, None] + env_means[None, :] - grand_mean
    residuals = Y - predicted
    v_gxe = float(np.nanvar(residuals[observed]))
    # V_total
    v_total = float(np.nanvar(Y[observed]))

    # Heritability across environments
    h2_broad = v_g / v_total if v_total > 0 else 0.0
    # G×E ratio
    gxe_ratio = v_gxe / v_total if v_total > 0 else 0.0
    # Type B genetic correlation: rB = σ²_G / (σ²_G + σ²_GxE)
    r_b = v_g / (v_g + v_gxe) if (v_g + v_gxe) > 0 else 0.0

    # Per-environment heritability (variance of genotype effects within env)
    h2_per_env = {}
    for j, env_id in enumerate(environments):
        col = Y[:, j]
        valid = ~np.isnan(col)
        if np.sum(valid) > 5:
            h2_per_env[env_id] = float(np.var(col[valid]) / v_total) if v_total > 0 else 0.0
        else:
            h2_per_env[env_id] = None

    if verbose:
        print(f"=== G×E Decomposition ===")
        print(f"  Trial: {trial_id}, Trait: {trait}")
        print(f"  {G} genotypes × {E} environments, "
              f"{n_obs}/{n_total} cells observed ({balance_ratio:.1%})")
        print(f"\n  Variance components:")
        print(f"    V_G   = {v_g:.4f}  ({v_g/v_total*100:.1f}%)")
        print(f"    V_E   = {v_e:.4f}  ({v_e/v_total*100:.1f}%)")
        print(f"    V_GxE = {v_gxe:.4f}  ({v_gxe/v_total*100:.1f}%)")
        print(f"    V_tot = {v_total:.4f}")
        print(f"\n  h²_broad = {h2_broad:.4f}")
        print(f"  G×E ratio = {gxe_ratio:.4f}")
        print(f"  Type B r_G = {r_b:.4f}")

    return {
        "trial_id": trial_id,
        "trait": trait,
        "v_g": v_g, "v_e": v_e, "v_gxe": v_gxe, "v_total": v_total,
        "h2_broad": h2_broad,
        "gxe_ratio": gxe_ratio,
        "type_b_correlation": r_b,
        "h2_per_environment": h2_per_env,
        "n_genotypes": G, "n_environments": E,
        "n_observed": n_obs, "n_missing": n_total - n_obs,
        "balance_ratio": balance_ratio,
    }


# ===================================================================
# M7. Design Diagnostics
# ===================================================================

def design_diagnostics(conn: GraphGWASConnection,
                       trial_id: str,
                       verbose: bool = True) -> dict:
    """Graph connectivity analysis of the experimental design.

    Builds bipartite graph (genotypes ↔ environments) and computes:
    - Balance ratio (observed / total cells)
    - Connected components (disconnected subgraphs = confounding)
    - Bridge genotypes (high betweenness = reference set)
    - Algebraic connectivity (Fiedler value = design informativeness)
    - Per-genotype environmental coverage
    - Per-environment genetic diversity

    Args:
        conn: database connection.
        trial_id: trial identifier.
        verbose: print progress.

    Returns:
        Comprehensive design diagnostics dict.
    """
    # Query the bipartite structure
    result = conn.execute_read(
        """
        MATCH (g:Genotype)-[:OBSERVED_IN]->(e:Environment {trial_id: $tid})
        RETURN g.id AS geno, e.id AS env
        """,
        {"tid": trial_id},
    )

    edges = [(r["geno"], r["env"]) for r in result]
    if not edges:
        return {"error": "No observations found for this trial"}

    genotypes = sorted(set(g for g, _ in edges))
    environments = sorted(set(e for _, e in edges))
    G, E = len(genotypes), len(environments)

    # Build bipartite graph
    B = nx.Graph()
    B.add_nodes_from(genotypes, bipartite="genotype")
    B.add_nodes_from(environments, bipartite="environment")
    B.add_edges_from(edges)

    n_obs = len(set(edges))  # unique genotype-environment combinations
    n_total = G * E
    balance_ratio = n_obs / n_total if n_total > 0 else 0

    # Connected components
    components = list(nx.connected_components(B))
    n_components = len(components)

    # Per-genotype coverage
    geno_coverage = {}
    for g in genotypes:
        n_envs = len(set(e for gg, e in edges if gg == g))
        geno_coverage[g] = n_envs

    # Per-environment diversity
    env_diversity = {}
    for e in environments:
        n_genos = len(set(g for g, ee in edges if ee == e))
        env_diversity[e] = n_genos

    # Bridge genotypes (high betweenness centrality)
    if B.number_of_edges() > 0 and n_components == 1:
        betweenness = nx.betweenness_centrality(B)
        # Filter to genotype nodes only
        geno_betweenness = {n: v for n, v in betweenness.items()
                           if n in genotypes}
        bridge_genotypes = sorted(geno_betweenness.items(),
                                  key=lambda x: x[1], reverse=True)[:20]
    else:
        bridge_genotypes = []

    # Algebraic connectivity (on genotype projection)
    # Project bipartite to genotype-genotype (shared environments)
    if G > 1 and E > 1:
        try:
            geno_proj = nx.bipartite.projected_graph(B, genotypes)
            if geno_proj.number_of_edges() > 0:
                L = nx.laplacian_matrix(geno_proj).toarray().astype(float)
                eigenvalues = np.sort(np.linalg.eigvalsh(L))
                fiedler = float(eigenvalues[1]) if len(eigenvalues) > 1 else 0.0
            else:
                fiedler = 0.0
        except Exception:
            fiedler = 0.0
    else:
        fiedler = 0.0

    if verbose:
        print(f"=== Design Diagnostics ===")
        print(f"  Trial: {trial_id}")
        print(f"  {G} genotypes × {E} environments")
        print(f"  Observations: {n_obs} / {n_total} ({balance_ratio:.1%} balanced)")
        print(f"  Connected components: {n_components}")
        print(f"  Algebraic connectivity (Fiedler): {fiedler:.4f}")
        print(f"  Genotype coverage: mean={np.mean(list(geno_coverage.values())):.1f} envs, "
              f"range=[{min(geno_coverage.values())}, {max(geno_coverage.values())}]")
        print(f"  Environment diversity: mean={np.mean(list(env_diversity.values())):.0f} genotypes, "
              f"range=[{min(env_diversity.values())}, {max(env_diversity.values())}]")
        if bridge_genotypes:
            print(f"  Top bridge genotypes (reference set candidates):")
            for g, bc in bridge_genotypes[:5]:
                print(f"    {g}: betweenness={bc:.4f}, "
                      f"coverage={geno_coverage.get(g, 0)} envs")

    return {
        "trial_id": trial_id,
        "n_genotypes": G,
        "n_environments": E,
        "n_observations": n_obs,
        "n_missing": n_total - n_obs,
        "balance_ratio": balance_ratio,
        "n_connected_components": n_components,
        "algebraic_connectivity": fiedler,
        "bridge_genotypes": [{"genotype": g, "betweenness": bc}
                             for g, bc in bridge_genotypes],
        "genotype_coverage": geno_coverage,
        "environment_diversity": env_diversity,
        "mean_genotype_coverage": float(np.mean(list(geno_coverage.values()))),
        "mean_env_diversity": float(np.mean(list(env_diversity.values()))),
    }


# ===================================================================
# M8. Graph-Based Imputation of Missing G×E Cells
# ===================================================================

def graph_impute_missing(conn: GraphGWASConnection,
                         trial_id: str,
                         trait: str,
                         genotype_id: str,
                         environment_id: str,
                         k_neighbors: int = 10,
                         verbose: bool = True) -> dict:
    """Predict missing genotype-environment phenotype via graph neighbors.

    Finds the k genotypes most similar to the target (by shared environment
    performance pattern) that WERE tested in the target environment.
    Weighted average of their trait values = imputed value.

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name.
        genotype_id: genotype to impute for.
        environment_id: environment to impute in.
        k_neighbors: number of nearest neighbors.
        verbose: print progress.

    Returns:
        dict with imputed_value, confidence, and neighbor details.
    """
    safe_trait = f"pheno_{trait.replace(' ', '_').replace('-', '_')}"

    # Get target genotype's performance in other environments
    target_result = conn.execute_read(
        f"""
        MATCH (g:Genotype {{id: $gid}})-[r:OBSERVED_IN]->(e:Environment {{trial_id: $tid}})
        WHERE r.{safe_trait} IS NOT NULL
        RETURN e.id AS env, r.{safe_trait} AS val
        """,
        {"gid": genotype_id, "tid": trial_id},
    )
    target_profile = {r["env"]: r["val"] for r in target_result}

    if not target_profile:
        return {"error": f"No observations for genotype {genotype_id}"}

    # Get all genotypes tested in the target environment
    candidate_result = conn.execute_read(
        f"""
        MATCH (g:Genotype)-[r:OBSERVED_IN]->(e:Environment {{id: $eid}})
        WHERE r.{safe_trait} IS NOT NULL AND g.id <> $gid
        RETURN g.id AS geno, r.{safe_trait} AS target_val
        """,
        {"eid": environment_id, "gid": genotype_id},
    )

    candidates = {r["geno"]: r["target_val"] for r in candidate_result}
    if not candidates:
        return {"error": f"No genotypes found in environment {environment_id}"}

    # For each candidate, compute similarity based on shared env performance
    neighbors = []
    for cand_id, cand_target_val in candidates.items():
        cand_result = conn.execute_read(
            f"""
            MATCH (g:Genotype {{id: $gid}})-[r:OBSERVED_IN]->(e:Environment {{trial_id: $tid}})
            WHERE r.{safe_trait} IS NOT NULL
            RETURN e.id AS env, r.{safe_trait} AS val
            """,
            {"gid": cand_id, "tid": trial_id},
        )
        cand_profile = {r["env"]: r["val"] for r in cand_result}

        # Shared environments (excluding target)
        shared_envs = set(target_profile.keys()) & set(cand_profile.keys())
        shared_envs.discard(environment_id)

        if len(shared_envs) < 2:
            continue

        t_vals = np.array([target_profile[e] for e in shared_envs])
        c_vals = np.array([cand_profile[e] for e in shared_envs])

        if np.std(t_vals) == 0 or np.std(c_vals) == 0:
            continue

        similarity = float(np.corrcoef(t_vals, c_vals)[0, 1])
        if np.isnan(similarity):
            continue

        neighbors.append({
            "genotype": cand_id,
            "similarity": similarity,
            "observed_value": cand_target_val,
            "n_shared_envs": len(shared_envs),
        })

    # Sort by similarity, take top k
    neighbors.sort(key=lambda n: n["similarity"], reverse=True)
    neighbors = neighbors[:k_neighbors]

    if not neighbors:
        return {"error": "No suitable neighbors found"}

    # Weighted imputation (similarity-weighted average)
    weights = np.array([max(n["similarity"], 0) for n in neighbors])
    values = np.array([n["observed_value"] for n in neighbors])

    if np.sum(weights) > 0:
        imputed = float(np.average(values, weights=weights))
    else:
        imputed = float(np.mean(values))

    # Confidence: mean similarity × sqrt(n_neighbors) / max_possible
    mean_sim = float(np.mean([n["similarity"] for n in neighbors]))
    confidence = float(min(1.0, mean_sim * np.sqrt(len(neighbors)) / np.sqrt(k_neighbors)))

    if verbose:
        print(f"=== Graph Imputation ===")
        print(f"  Genotype: {genotype_id}, Environment: {environment_id}")
        print(f"  Imputed value: {imputed:.4f}")
        print(f"  Confidence: {confidence:.4f}")
        print(f"  Neighbors: {len(neighbors)}")
        for n in neighbors[:5]:
            print(f"    {n['genotype']}: sim={n['similarity']:.4f}, "
                  f"val={n['observed_value']:.4f}")

    return {
        "genotype": genotype_id,
        "environment": environment_id,
        "trait": trait,
        "imputed_value": imputed,
        "confidence": confidence,
        "n_neighbors": len(neighbors),
        "neighbors": neighbors,
    }


# ===================================================================
# M9. Integrated MET Report
# ===================================================================

def met_report(conn: GraphGWASConnection,
               trial_id: str,
               trait: str,
               verbose: bool = True) -> dict:
    """Run comprehensive MET analysis: diagnostics + similarity + G×E.

    Args:
        conn: database connection.
        trial_id: trial identifier.
        trait: trait name.
        verbose: print progress.

    Returns:
        Unified report with all MET analysis results.
    """
    if verbose:
        print("=" * 60)
        print(f"  MET Report: {trial_id} / {trait}")
        print("=" * 60)

    results = {}

    # Design diagnostics
    if verbose:
        print("\n" + "-" * 60)
    results["diagnostics"] = design_diagnostics(conn, trial_id, verbose=verbose)

    # G×E decomposition
    if verbose:
        print("\n" + "-" * 60)
    results["gxe"] = gxe_decomposition(conn, trial_id, trait, verbose=verbose)

    # Environment similarity + mega-environments
    if verbose:
        print("\n" + "-" * 60)
    results["mega_environments"] = mega_environments(
        conn, trial_id, trait, verbose=verbose,
    )

    # Summary
    if verbose:
        diag = results["diagnostics"]
        gxe = results["gxe"]
        mega = results["mega_environments"]
        print("\n" + "=" * 60)
        print("  Summary")
        print("=" * 60)
        print(f"  Design: {diag['n_genotypes']} genotypes × "
              f"{diag['n_environments']} environments "
              f"({diag['balance_ratio']:.1%} balanced)")
        print(f"  Components: {diag['n_connected_components']}, "
              f"Fiedler: {diag['algebraic_connectivity']:.4f}")
        if "v_g" in gxe:
            print(f"  V_G={gxe['v_g']:.4f}, V_E={gxe['v_e']:.4f}, "
                  f"V_GxE={gxe['v_gxe']:.4f}")
            print(f"  h²_broad={gxe['h2_broad']:.4f}, "
                  f"Type B r_G={gxe['type_b_correlation']:.4f}")
        print(f"  Mega-environments: {mega['n_clusters']}")

    return results
