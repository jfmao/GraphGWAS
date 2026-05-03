"""Build a human catalogue extension from the Billmann et al. 2026 Cell
HAP1 genetic interaction map.

Strategy:
  - 71 nonredundant CORUM protein complexes (File_S11) define the
    pair set: each complex's first 2 genes (alphabetical) become the
    "anchor" pair representing that complex.
  - File_S22 "Network evaluations_CORUM" gives integrated/depmap/GI
    AUPRC per complex — used to tier each entry.
  - Everything goes under id "human.billmann.<complex_short_name>".

Schema matches the existing docs/groundtruth/epistasis_pairs.json human
arm (gene_1.hgnc, gene_1.locus, validation list, mechanism, etc.). The
locus field is left null for now — chromosomal coordinates will be
filled in when the human cache is wired (paper-#2-revision target;
see results_summary_y8 cache-remediation discussion).

Usage: python scripts/curate_billmann_human_pairs.py

Writes:
  - docs/groundtruth/billmann_human_pairs_curated.json (intermediate)
  - then updates docs/groundtruth/epistasis_pairs.json (appends to
    species.human.pairs after re-tiering existing 8 entries)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from collections import OrderedDict

import openpyxl

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data" / "billmann2026"
FILE_S11 = DATA_DIR / "File_S11.xlsx"
FILE_S22 = DATA_DIR / "File_S22.xlsx"
CATALOGUE = REPO / "docs" / "groundtruth" / "epistasis_pairs.json"
INTERMEDIATE = REPO / "docs" / "groundtruth" / "billmann_human_pairs_curated.json"

# ---------------------------------------------------------------------------
# Auto-download from the Boone lab supplement if files are missing
# ---------------------------------------------------------------------------

BOONE_BASE = "https://boonelab.ccbr.utoronto.ca/supplement/billmanncostanzo2026/assets/data_files"
NEEDED_FILES = {
    "File_S11.xlsx": "CORUM nonredundant complexes + module annotations",
    "File_S22.xlsx": "HAP1+DepMap integrated network evaluations (AUPRC per CORUM)",
}


def ensure_files() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    import urllib.request
    for fname, desc in NEEDED_FILES.items():
        path = DATA_DIR / fname
        if path.exists() and path.stat().st_size > 1000:
            continue
        url = f"{BOONE_BASE}/{fname}"
        print(f"Downloading {fname} ({desc}) from {url}...")
        urllib.request.urlretrieve(url, path)
        print(f"  → {path.stat().st_size/1e6:.2f} MB")


ensure_files()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def short_id(complex_name: str) -> str:
    """Make a stable ID slug from the complex name."""
    s = re.sub(r"\s+", "_", complex_name.strip())
    s = re.sub(r"[^A-Za-z0-9_-]", "", s)
    return s[:60]


def safe_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Load File_S11 nonredundant CORUM complexes
# ---------------------------------------------------------------------------

print(f"Reading {FILE_S11.name}...")
wb11 = openpyxl.load_workbook(FILE_S11, read_only=True, data_only=True)
ws = wb11["CORUM_nonredundant_complexes"]
nonred: list = []
for r in list(ws.iter_rows(values_only=True))[1:]:
    if r[1] and r[2]:
        genes = sorted({g.strip() for g in r[2].split(";") if g.strip()})
        if len(genes) < 2:
            continue
        nonred.append({
            "bioprocess_region": (r[0] or "").strip(),
            "complex_name": r[1].strip(),
            "genes": genes,
        })
wb11.close()
print(f"  {len(nonred)} nonredundant CORUM complexes (≥2 members)")

# ---------------------------------------------------------------------------
# Load File_S22 AUPRC evaluations + index by lowercase complex name
# ---------------------------------------------------------------------------

print(f"Reading {FILE_S22.name} → CORUM AUPRC evaluations...")
wb22 = openpyxl.load_workbook(FILE_S22, read_only=True, data_only=True)
ws = wb22["Network evaluations_CORUM"]
auprc_idx: dict = {}
for r in list(ws.iter_rows(values_only=True))[1:]:
    name = r[1]
    if not name:
        continue
    a = safe_float(r[5])
    if a is None:
        continue
    key = name.strip().lower()
    # Use the BEST (max) integrated AUPRC if multiple rows match same name
    prev = auprc_idx.get(key)
    if prev is None or a > prev["auprc_integrated"]:
        auprc_idx[key] = {
            "n_genes_eval": r[2],
            "auprc_depmap": safe_float(r[3]),
            "auprc_gi": safe_float(r[4]),
            "auprc_integrated": a,
        }
wb22.close()
print(f"  {len(auprc_idx)} unique CORUM names with numeric integrated AUPRC")

# Match nonredundant complexes
matched = 0
for c in nonred:
    a = auprc_idx.get(c["complex_name"].lower())
    if a:
        c.update(a)
        matched += 1
print(f"  matched {matched}/{len(nonred)} nonredundant complexes")

# ---------------------------------------------------------------------------
# Build catalogue entries
# ---------------------------------------------------------------------------

# Tier rule: integrated AUPRC ≥ 0.5 OR (gi_auprc ≥ 0.4 AND depmap_auprc ≥ 0.4)
# → "A_billmann_corum_high"; else "B_billmann_corum_supported";
# unmatched → "C_billmann_corum_only" (Billmann GI evidence only, no DepMap eval).
def tier(c):
    a = c.get("auprc_integrated")
    if a is None:
        return "C_billmann_corum_only"
    if a >= 0.5:
        return "A_billmann_corum_high"
    if (c.get("auprc_gi") or 0) >= 0.4 and (c.get("auprc_depmap") or 0) >= 0.4:
        return "A_billmann_corum_high"
    return "B_billmann_corum_supported"


# Anchor: pick the FA core complex (Fanconi anemia) as the human anchor
# for this tier — clinically famous AND in the AUPRC top-15.
ANCHOR_COMPLEX_NAMES = {
    "fa core complex (fanconi anemia core complex)",
    "fanconi anemia, core complex",
}

new_entries: list = []
for c in nonred:
    g1, g2 = c["genes"][0], c["genes"][1]  # first 2 alphabetical
    cid = short_id(c["complex_name"])
    is_anchor = c["complex_name"].lower() in ANCHOR_COMPLEX_NAMES
    bp = c.get("bioprocess_region") or "unannotated"

    validation_strs = [
        "HAP1 combinatorial CRISPR-CRISPR screen "
        "(Billmann et al. 2026 Cell, ~4M gene-pairs tested, qGI+FDR)",
        f"CORUM nonredundant protein complex co-membership: \"{c['complex_name']}\"",
    ]
    if "auprc_integrated" in c:
        dm = c.get("auprc_depmap"); gi = c.get("auprc_gi")
        dm_s = f"{dm:.3f}" if dm is not None else "n/a"
        gi_s = f"{gi:.3f}" if gi is not None else "n/a"
        validation_strs.append(
            f"HAP1+DepMap integrated AUPRC = {c['auprc_integrated']:.3f} "
            f"(DepMap-only: {dm_s}; GI-only: {gi_s})"
        )

    entry = OrderedDict([
        ("id", f"human.billmann.{cid}"),
        ("gene_1", {"hgnc": g1, "locus": None}),
        ("gene_2", {"hgnc": g2, "locus": None}),
        ("phenotype",
         f"Cellular fitness in HAP1 (proxy for cell viability/proliferation)"),
        ("validation", validation_strs),
        ("panel_natural_variation", True),
        ("anchor", is_anchor),
        ("mechanism",
         f"Co-members of CORUM complex \"{c['complex_name']}\" "
         f"(bioprocess: {bp}); double-perturbation in HAP1 produces a "
         f"detectable qGI score, integrated with DepMap "
         f"co-essentiality where AUPRC available."),
        ("primary_citation",
         "Billmann et al. 2026 Cell doi:10.1016/j.cell.2026.03.044; "
         "Tsitsiridis et al. 2023 NAR doi:10.1093/nar/gkac1015 (CORUM 4.0)"),
        ("expected_motif", ["protein_interaction"]),
        ("tier", tier(c)),
        ("billmann_evidence", {
            "complex_name": c["complex_name"],
            "bioprocess_region": bp,
            "n_complex_members": len(c["genes"]),
            "all_complex_members": c["genes"],
            "auprc_integrated": c.get("auprc_integrated"),
            "auprc_depmap": c.get("auprc_depmap"),
            "auprc_gi": c.get("auprc_gi"),
        }),
    ])
    new_entries.append(entry)

# Save intermediate
INTERMEDIATE.write_text(json.dumps(new_entries, indent=1))
print(f"\nWrote {len(new_entries)} curated entries → {INTERMEDIATE}")

# Tier distribution
from collections import Counter
tdist = Counter(e["tier"] for e in new_entries)
for t, n in sorted(tdist.items()):
    print(f"  tier {t}: {n}")

# ---------------------------------------------------------------------------
# Update epistasis_pairs.json — append to species.human.pairs
# ---------------------------------------------------------------------------

print(f"\nUpdating {CATALOGUE}...")
catalogue = json.loads(CATALOGUE.read_text())
human_section = catalogue["species"]["human"]
existing = human_section["pairs"]
print(f"  existing human pairs: {len(existing)}")
# Re-tier existing 8 entries to "A_pharmacogenetic_or_disease_modifier"
# so the two tier families are clearly distinct
for p in existing:
    old_tier = p.get("tier", "?")
    if old_tier in {"A", "B", "C"}:
        p["tier"] = "A_gwas_pharmacogenetic" if old_tier == "A" else "B_disease_modifier"

human_section["pairs"] = existing + new_entries
human_section["counts"] = {
    "n_total": len(human_section["pairs"]),
    "n_gwas_pharmacogenetic": len(existing),
    "n_billmann_corum": len(new_entries),
}
human_section["tier_legend"] = {
    "A_gwas_pharmacogenetic":
        "Population-scale GWAS-validated drug-response or disease-modifier "
        "epistasis pair (HLA-B×ERAP1, TPMT×NUDT15, VKORC1×CYP2C9, etc.)",
    "B_disease_modifier":
        "Disease-modifier pair with mechanistic evidence but limited "
        "population-scale replication.",
    "A_billmann_corum_high":
        "CORUM protein complex co-members; HAP1 GI screen + DepMap "
        "co-essentiality both confirm functional grouping (integrated "
        "AUPRC ≥ 0.5 OR both component AUPRCs ≥ 0.4).",
    "B_billmann_corum_supported":
        "CORUM protein complex co-members; HAP1 GI screen + DepMap "
        "evaluated but integrated AUPRC < 0.5.",
    "C_billmann_corum_only":
        "CORUM protein complex co-members; not in DepMap AUPRC eval set "
        "(Billmann HAP1 GI evidence only).",
}

CATALOGUE.write_text(json.dumps(catalogue, indent=1))
print(f"  human pairs after merge: {len(human_section['pairs'])}")
print(f"  Wrote {CATALOGUE}")
