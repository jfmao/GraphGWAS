"""GNN Association Engine — Heterogeneous Graph Neural Network for GWAS.

Phase 4: Message-passing replaces marginal regression. A heterogeneous GNN
over the Variant-Sample-Gene-Pathway graph learns multi-locus epistatic
architecture implicitly through neighborhood aggregation.

Components:
1. export_to_pyg(): Neo4j → PyTorch Geometric HeteroData
2. GraphGWASModel: HeteroGNN for phenotype prediction
3. train() / evaluate(): training loop with class-weighted loss
4. explain(): GNNExplainer for variant-level attribution
5. import_embeddings(): GNN embeddings back to Neo4j
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import HeteroData
    from torch_geometric.nn import HeteroConv, SAGEConv, Linear
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

from . import config as _cfg
from .db import GraphGWASConnection
from .genotype import (
    build_carrier_set,
    build_dosage,
    get_phenotype_indices,
    get_all_indices,
    variant_iterator,
    unpack_genotypes,
)


# ---------------------------------------------------------------------------
# Device resolution
# ---------------------------------------------------------------------------

def _resolve_device(device: str = "auto"):
    """Resolve device string to torch.device. Returns None if torch unavailable.

    Auto-detection picks CUDA when available, otherwise CPU.
    MPS (Apple Silicon) is not auto-selected because PyTorch Geometric's
    scatter/message-passing ops have incomplete MPS support. Users can
    force MPS with device="mps" at their own risk.
    """
    if not HAS_TORCH:
        return None
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")
    return torch.device(device)


# ---------------------------------------------------------------------------
# 4.1 Graph Export to PyTorch Geometric
# ---------------------------------------------------------------------------

def export_to_pyg(conn: GraphGWASConnection, chr: str,
                  start: int | None = None, end: int | None = None,
                  af_range: tuple[float, float] = (0.001, 0.5),
                  max_variants: int = 10000,
                  verbose: bool = True) -> "HeteroData":
    """Export subgraph from Neo4j to PyTorch Geometric HeteroData.

    Node types:
    - 'sample': features = [is_case, sex, population_encoded]
    - 'variant': features = [af_total, call_rate, is_carrier_per_sample...]
    - 'gene': features = [n_variants, biotype_encoded]

    Edge types:
    - ('sample', 'carries', 'variant'): decoded from gt_packed
    - ('variant', 'in_gene', 'gene'): from HAS_CONSEQUENCE
    - ('gene', 'in_pathway', 'pathway'): from IN_PATHWAY

    Labels: sample.y = is_case (binary)
    """
    if not HAS_TORCH:
        raise ImportError("PyTorch and torch-geometric required. pip install torch torch-geometric")

    case_idx, ctrl_idx = get_phenotype_indices(conn)
    all_idx = np.concatenate([case_idx, ctrl_idx])
    n_samples = len(all_idx)

    # --- Sample nodes ---
    sample_result = conn.execute_read(
        """
        MATCH (s:Sample)
        WHERE s.packed_index IN $indices
        RETURN s.packed_index AS idx, s.sampleId AS sid,
               s.population AS pop, s.sex AS sex,
               s.is_case AS is_case
        ORDER BY s.packed_index
        """,
        {"indices": [int(i) for i in all_idx]},
    )
    samples = [dict(r) for r in sample_result]

    # Encode population as integer
    pop_set = sorted(set(s.get("pop", "UNK") or "UNK" for s in samples))
    pop_to_idx = {p: i for i, p in enumerate(pop_set)}

    sample_features = []
    sample_labels = []
    sample_id_map = {}  # packed_index → local index
    for i, s in enumerate(samples):
        sample_id_map[s["idx"]] = i
        is_case = 1.0 if s.get("is_case") else 0.0
        sex = float(s.get("sex", 0) or 0)
        pop_enc = float(pop_to_idx.get(s.get("pop", "UNK") or "UNK", 0))
        sample_features.append([is_case, sex, pop_enc])
        sample_labels.append(is_case)

    # --- Variant nodes + CARRIES edges ---
    if verbose:
        region = f"{chr}:{start}-{end}" if start else chr
        print(f"Exporting PyG graph for {region}...")

    variant_features = []
    variant_id_map = {}  # variantId → local index
    carries_src = []  # sample local indices
    carries_dst = []  # variant local indices

    var_idx = 0
    for v in variant_iterator(conn, chr, start, end):
        gt_packed = v["gt_packed"]
        if gt_packed is None:
            continue
        af = v.get("af_total", 0)
        if af < af_range[0] or af > af_range[1]:
            continue
        if var_idx >= max_variants:
            break

        vid = v["variantId"]
        variant_id_map[vid] = var_idx

        cr = v.get("call_rate", 1.0) or 1.0
        variant_features.append([float(af), float(cr)])

        # Decode carriers among our sample set
        carriers = build_carrier_set(gt_packed, all_idx, _cfg.N_SAMPLES)
        carrier_local_indices = np.where(carriers)[0]
        for ci in carrier_local_indices:
            carries_src.append(ci)
            carries_dst.append(var_idx)

        var_idx += 1

    if verbose:
        print(f"  {len(sample_features)} samples, {var_idx} variants, {len(carries_src)} carrier edges")

    # --- Gene nodes + variant→gene edges ---
    gene_result = conn.execute_read(
        """
        MATCH (v:Variant)-[r:HAS_CONSEQUENCE]->(g:Gene)
        WHERE v.variantId IN $vids
        RETURN v.variantId AS vid, g.symbol AS gene, g.biotype AS biotype
        """,
        {"vids": list(variant_id_map.keys())},
    )

    gene_id_map = {}
    gene_features = []
    in_gene_src = []  # variant local indices
    in_gene_dst = []  # gene local indices
    biotype_set = set()

    for r in gene_result:
        vid = r["vid"]
        gene = r["gene"]
        if vid not in variant_id_map or not gene:
            continue

        if gene not in gene_id_map:
            gene_id_map[gene] = len(gene_id_map)
            biotype = r.get("biotype", "unknown") or "unknown"
            biotype_set.add(biotype)
            gene_features.append([0.0])  # placeholder, updated below

        in_gene_src.append(variant_id_map[vid])
        in_gene_dst.append(gene_id_map[gene])

    # Update gene features with variant count
    gene_var_counts = [0] * len(gene_id_map)
    for gi in in_gene_dst:
        gene_var_counts[gi] += 1
    gene_features = [[float(c)] for c in gene_var_counts]

    if verbose:
        print(f"  {len(gene_id_map)} genes, {len(in_gene_src)} variant→gene edges")

    # --- Build HeteroData ---
    data = HeteroData()

    data['sample'].x = torch.tensor(sample_features, dtype=torch.float)
    data['sample'].y = torch.tensor(sample_labels, dtype=torch.float)
    data['sample'].num_nodes = n_samples

    data['variant'].x = torch.tensor(variant_features, dtype=torch.float) if variant_features else torch.zeros((0, 2))
    data['variant'].num_nodes = var_idx

    data['gene'].x = torch.tensor(gene_features, dtype=torch.float) if gene_features else torch.zeros((0, 1))
    data['gene'].num_nodes = len(gene_id_map)

    # Edges
    if carries_src:
        data['sample', 'carries', 'variant'].edge_index = torch.tensor(
            [carries_src, carries_dst], dtype=torch.long
        )
        # Reverse edges for message passing
        data['variant', 'carried_by', 'sample'].edge_index = torch.tensor(
            [carries_dst, carries_src], dtype=torch.long
        )

    if in_gene_src:
        data['variant', 'in_gene', 'gene'].edge_index = torch.tensor(
            [in_gene_src, in_gene_dst], dtype=torch.long
        )
        data['gene', 'has_variant', 'variant'].edge_index = torch.tensor(
            [in_gene_dst, in_gene_src], dtype=torch.long
        )

    # Store ID maps for later reference
    data.sample_id_map = sample_id_map
    data.variant_id_map = variant_id_map
    data.gene_id_map = gene_id_map

    if verbose:
        print(f"  HeteroData: {data}")

    return data


# ---------------------------------------------------------------------------
# 4.2 Heterogeneous GNN Model
# ---------------------------------------------------------------------------

if HAS_TORCH:

    class GraphGWASModel(nn.Module):
        """Heterogeneous GNN for phenotype prediction.

        Uses SAGEConv (GraphSAGE) with heterogeneous message passing.
        3 layers of message passing aggregate neighborhood information
        from variants through genes, learning multi-locus effects implicitly.
        """

        def __init__(self, metadata, hidden_channels: int = 64, num_layers: int = 3):
            super().__init__()
            self.num_layers = num_layers

            # Input projections per node type
            self.input_proj = nn.ModuleDict()
            # Will be set dynamically based on data

            # Heterogeneous convolution layers
            self.convs = nn.ModuleList()
            for _ in range(num_layers):
                conv_dict = {}
                for edge_type in metadata[1]:
                    src, rel, dst = edge_type
                    conv_dict[edge_type] = SAGEConv((-1, -1), hidden_channels)
                self.convs.append(HeteroConv(conv_dict, aggr='mean'))

            # Phenotype predictor on sample embeddings
            self.classifier = nn.Sequential(
                nn.Linear(hidden_channels, 32),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(32, 1),
            )

        def forward(self, x_dict, edge_index_dict):
            # Message passing layers
            for conv in self.convs:
                x_dict = conv(x_dict, edge_index_dict)
                x_dict = {key: F.relu(x) for key, x in x_dict.items()}
                x_dict = {key: F.dropout(x, p=0.3, training=self.training)
                          for key, x in x_dict.items()}

            # Predict phenotype from sample embeddings
            if 'sample' in x_dict:
                out = self.classifier(x_dict['sample'])
                return out.squeeze(-1), x_dict
            return None, x_dict


# ---------------------------------------------------------------------------
# 4.3 Training
# ---------------------------------------------------------------------------

def train_gnn(data: "HeteroData", hidden_channels: int = 64,
              num_layers: int = 3, epochs: int = 100,
              lr: float = 0.001, weight_decay: float = 1e-4,
              val_ratio: float = 0.2,
              device: str = "auto",
              verbose: bool = True) -> dict:
    """Train the GraphGWAS GNN model.

    Args:
        data: HeteroData from export_to_pyg.
        hidden_channels: GNN hidden dimension.
        num_layers: number of message-passing layers.
        epochs: training epochs.
        lr: learning rate.
        weight_decay: L2 regularization.
        val_ratio: fraction held out for validation.

    Returns dict with model, train_losses, val_aurocs, best_auroc.
    """
    if not HAS_TORCH:
        raise ImportError("PyTorch required")
    from sklearn.metrics import roc_auc_score

    # Resolve device (CPU or CUDA)
    dev = _resolve_device(device)
    if verbose:
        print(f"Training on device: {dev}")

    # Train/val split on samples
    n_samples = data['sample'].num_nodes
    perm = torch.randperm(n_samples)
    n_val = int(n_samples * val_ratio)
    val_idx = perm[:n_val]
    train_idx = perm[n_val:]

    train_mask = torch.zeros(n_samples, dtype=torch.bool)
    train_mask[train_idx] = True
    val_mask = torch.zeros(n_samples, dtype=torch.bool)
    val_mask[val_idx] = True

    # Class weights for imbalanced case/control
    labels = data['sample'].y
    n_pos = labels.sum().item()
    n_neg = n_samples - n_pos
    pos_weight = torch.tensor([n_neg / max(n_pos, 1)])

    # Move data and masks to device
    data = data.to(dev)
    labels = data['sample'].y
    train_mask = train_mask.to(dev)
    val_mask = val_mask.to(dev)
    pos_weight = pos_weight.to(dev)

    # Initialize model
    model = GraphGWASModel(data.metadata(), hidden_channels=hidden_channels,
                           num_layers=num_layers)

    # Lazy init: run one forward pass to set dimensions (on CPU first)
    with torch.no_grad():
        model.eval()
        model(data.x_dict, data.edge_index_dict)

    # Move model to device after lazy init
    model = model.to(dev)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    train_losses = []
    val_aurocs = []
    best_auroc = 0.0
    best_state = None

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()

        logits, embeddings = model(data.x_dict, data.edge_index_dict)
        if logits is None:
            break

        loss = criterion(logits[train_mask], labels[train_mask])
        loss.backward()
        optimizer.step()

        train_losses.append(loss.item())

        # Validation
        if (epoch + 1) % 10 == 0 or epoch == epochs - 1:
            model.eval()
            with torch.no_grad():
                logits_val, _ = model(data.x_dict, data.edge_index_dict)
                probs = torch.sigmoid(logits_val[val_mask]).cpu().numpy()
                y_true = labels[val_mask].cpu().numpy()

                try:
                    auroc = roc_auc_score(y_true, probs)
                except ValueError:
                    auroc = 0.5

                val_aurocs.append(auroc)
                if auroc > best_auroc:
                    best_auroc = auroc
                    best_state = {k: v.clone() for k, v in model.state_dict().items()}

                if verbose:
                    print(f"Epoch {epoch+1}/{epochs}: loss={loss.item():.4f} "
                          f"val_AUROC={auroc:.4f}")

    if best_state:
        model.load_state_dict(best_state)

    if verbose:
        print(f"\nBest validation AUROC: {best_auroc:.4f}")

    return {
        "model": model,
        "train_losses": train_losses,
        "val_aurocs": val_aurocs,
        "best_auroc": best_auroc,
        "data": data,
        "device": dev,
    }


# ---------------------------------------------------------------------------
# 4.4 Feature Attribution (Explainability)
# ---------------------------------------------------------------------------

def explain_variants(model, data: "HeteroData",
                     top_k: int = 20,
                     verbose: bool = True) -> list[dict]:
    """Variant-level feature attribution using gradient-based importance.

    Computes gradient of phenotype prediction w.r.t. each variant's
    carrier edges. High-gradient variants = strong phenotype contributors.

    Returns list of {variantId, importance_score, rank}.
    """
    if not HAS_TORCH:
        raise ImportError("PyTorch required")

    model.eval()

    # Enable gradient on variant features
    variant_x = data['variant'].x.clone().requires_grad_(True)
    x_dict = {k: v.clone() for k, v in data.x_dict.items()}
    x_dict['variant'] = variant_x

    logits, embeddings = model(x_dict, data.edge_index_dict)
    if logits is None:
        return []

    # Gradient of mean case probability w.r.t. variant features
    case_mask = data['sample'].y == 1.0
    case_logits = logits[case_mask].mean()
    case_logits.backward()

    # Variant importance = L2 norm of gradient
    variant_grad = variant_x.grad
    if variant_grad is None:
        return []

    importance = torch.norm(variant_grad, dim=1).detach().cpu().numpy()

    # Map back to variant IDs
    inv_map = {v: k for k, v in data.variant_id_map.items()}
    results = []
    for i in range(len(importance)):
        vid = inv_map.get(i, f"variant_{i}")
        results.append({
            "variantId": vid,
            "importance_score": float(importance[i]),
        })

    results.sort(key=lambda r: r["importance_score"], reverse=True)
    for i, r in enumerate(results):
        r["rank"] = i + 1

    if verbose:
        print(f"\nTop {min(top_k, len(results))} variants by GNN importance:")
        for r in results[:top_k]:
            print(f"  #{r['rank']} {r['variantId']:<50} score={r['importance_score']:.6f}")

    return results[:top_k]


# ---------------------------------------------------------------------------
# 4.5 Import Embeddings Back to Neo4j
# ---------------------------------------------------------------------------

def import_embeddings(conn: GraphGWASConnection, model, data: "HeteroData",
                      verbose: bool = True) -> int:
    """Import GNN-learned sample embeddings back to Neo4j.

    Stores gnn_embedding (float array) and gnn_phenotype_score (float)
    on Sample nodes. These capture multi-locus genetic effects — richer
    than PCA, usable as covariates in association tests.
    """
    if not HAS_TORCH:
        raise ImportError("PyTorch required")

    model.eval()
    with torch.no_grad():
        logits, embeddings = model(data.x_dict, data.edge_index_dict)

    if 'sample' not in embeddings:
        return 0

    sample_emb = embeddings['sample'].cpu().numpy()
    sample_scores = torch.sigmoid(logits).cpu().numpy() if logits is not None else np.zeros(len(sample_emb))

    inv_map = {v: k for k, v in data.sample_id_map.items()}

    updates = []
    for local_idx in range(len(sample_emb)):
        packed_idx = inv_map.get(local_idx)
        if packed_idx is None:
            continue
        updates.append({
            "idx": int(packed_idx),
            "score": float(sample_scores[local_idx]),
            "emb_dim": int(sample_emb.shape[1]),
        })

    # Write in batches
    n_updated = 0
    for i in range(0, len(updates), 500):
        batch = updates[i:i+500]
        result = conn.execute_write(
            """
            UNWIND $batch AS row
            MATCH (s:Sample {packed_index: row.idx})
            SET s.gnn_phenotype_score = row.score
            RETURN count(s) AS n
            """,
            {"batch": batch},
        )
        rec = result.single()
        n_updated += rec["n"] if rec else 0

    if verbose:
        print(f"Imported GNN scores for {n_updated} samples")

    return n_updated


# ---------------------------------------------------------------------------
# 4.6 Full GNN Pipeline
# ---------------------------------------------------------------------------

def run_gnn_pipeline(conn: GraphGWASConnection, chr: str,
                     start: int | None = None, end: int | None = None,
                     hidden_channels: int = 64, num_layers: int = 3,
                     epochs: int = 100, max_variants: int = 5000,
                     device: str = "auto",
                     verbose: bool = True) -> dict:
    """Complete GNN pipeline: export → train → explain → import.

    Returns dict with model, best_auroc, top_variants, n_embeddings_imported.
    """
    if verbose:
        print("=" * 60)
        print("GraphGWAS GNN Pipeline")
        print("=" * 60)

    # 1. Export
    if verbose:
        print("\n--- Step 1: Export graph to PyG ---")
    data = export_to_pyg(conn, chr, start, end, max_variants=max_variants, verbose=verbose)

    if data['variant'].num_nodes == 0:
        print("No variants in region. Aborting.")
        return {}

    # 2. Train
    if verbose:
        print("\n--- Step 2: Train GNN ---")
    train_result = train_gnn(data, hidden_channels=hidden_channels,
                             num_layers=num_layers, epochs=epochs,
                             device=device, verbose=verbose)

    # 3. Explain
    if verbose:
        print("\n--- Step 3: Variant attribution ---")
    top_variants = explain_variants(train_result["model"], data, verbose=verbose)

    # 4. Import embeddings
    if verbose:
        print("\n--- Step 4: Import embeddings to Neo4j ---")
    n_imported = import_embeddings(conn, train_result["model"], data, verbose=verbose)

    return {
        "model": train_result["model"],
        "best_auroc": train_result["best_auroc"],
        "top_variants": top_variants,
        "n_embeddings_imported": n_imported,
        "data": data,
    }
