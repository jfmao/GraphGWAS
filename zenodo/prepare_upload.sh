#!/usr/bin/env bash
# prepare_upload.sh — stage all GraphGWAS artifacts for Zenodo deposit
#
# Run this from the project root ONCE at paper-acceptance time:
#     cd /path/to/GraphGWAS && bash zenodo/prepare_upload.sh
#
# It populates every zenodo/ subfolder from live project data and writes
# MANIFEST.md with SHA-256 hashes for upload verification.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ZEN="$ROOT/zenodo"

echo "=== GraphGWAS Zenodo staging ==="
echo "Project root: $ROOT"
echo "Zenodo dir:   $ZEN"

# ===================================================================
# 1. Graph-database dumps
# ===================================================================
echo ""
echo "[1/6] Copying graph-database dumps..."
if [ -f "$ROOT/backups/yeast_1011_1.92M_variants_with_annotations.dump" ]; then
    cp "$ROOT/backups/yeast_1011_1.92M_variants_with_annotations.dump" \
       "$ZEN/graph_dumps/yeast_1011_v0.1.dump"
    echo "  ✓ yeast dump"
else
    echo "  ⚠ yeast dump not at expected path"
fi
if [ -f "$ROOT/backups/human_1kg_70.7M_multiomics_STRING_eQTL_conservation.dump" ]; then
    cp "$ROOT/backups/human_1kg_70.7M_multiomics_STRING_eQTL_conservation.dump" \
       "$ZEN/graph_dumps/human_1kg_multiomics_v0.1.dump"
    echo "  ✓ human dump"
else
    echo "  ⚠ human dump not at expected path"
fi

# ===================================================================
# 2. Benchmark JSON outputs
# ===================================================================
echo ""
echo "[2/6] Copying benchmark outputs..."
BENCH_SRC="$ROOT/results/benchmark_v2"
BENCH_DST="$ZEN/benchmark_outputs"
for f in \
    "weak_signal_l1_vs_susie_79rep" \
    "weak_signal_l1_vs_susie_100rep" \
    "wu2026_baselines_comparison" \
    "polyfun_proxy_comparison" \
    "pip_calibration_200rep" \
    "null_fpr_100rep" \
    "power_vs_N_30rep_per_N" \
    "cross_ancestry_1kg_30rep_per_ancestry"; do
    # Heuristic search — files may have slightly different paths
    matches=$(find "$BENCH_SRC" -maxdepth 3 -name "*${f//_/*}*.json" 2>/dev/null | head -3)
    if [ -n "$matches" ]; then
        for m in $matches; do
            cp "$m" "$BENCH_DST/" && echo "  ✓ $(basename "$m")"
        done
    else
        echo "  ⚠ missing: $f"
    fi
done

# ===================================================================
# 3. Pan-UKB per-locus results
# ===================================================================
echo ""
echo "[3/6] Copying Pan-UKB results..."
PANU_SRC="$ROOT/results/panukb"
if [ -d "$PANU_SRC" ]; then
    cp "$PANU_SRC"/*.json "$ZEN/panukb_results/" 2>/dev/null || true
    ls "$ZEN/panukb_results/" | grep -v README | head -5 | while read f; do
        echo "  ✓ $f"
    done
fi

# ===================================================================
# 4. Figure source data (extract per-figure CSVs)
# ===================================================================
echo ""
echo "[4/6] Extracting figure source data..."
for i in 1 2 3 4 5 6 7 8; do
    mkdir -p "$ZEN/figure_source_data/fig${i}"
done
echo "  (figure-specific CSV extraction to be run manually from the benchmark JSONs at this point)"

# ===================================================================
# 5. Simulation seeds
# ===================================================================
echo ""
echo "[5/6] Collating simulation seeds..."
SEED_DST="$ZEN/simulation_seeds"
python3 - <<'PYEOF'
import glob, json, os, csv
from pathlib import Path

ROOT = Path(os.environ.get("ROOT", "/mnt/data/GraphGWAS"))
DST = ROOT / "zenodo/simulation_seeds"
DST.mkdir(parents=True, exist_ok=True)

# Walk every benchmark JSON we can find and extract seed, window, causal.
for jf in glob.glob(str(ROOT / "results/benchmark_v2/**/*.json"), recursive=True):
    try:
        d = json.load(open(jf))
    except Exception:
        continue
    if not isinstance(d, dict) or "records" not in d:
        continue
    scenario = d.get("scenario", Path(jf).stem)
    tsv = DST / f"{scenario}.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["replicate_id", "seed", "chr", "start", "end", "causal_variant_id"])
        for rec in d.get("records", []):
            w.writerow([
                rec.get("replicate_id", ""),
                rec.get("seed", ""),
                rec.get("chr", ""),
                rec.get("start", ""),
                rec.get("end", ""),
                rec.get("causal_variant_id", ""),
            ])
    print(f"  ✓ {tsv.name}")
PYEOF

# ===================================================================
# 6a. Multi-omics graph caches (per species, baseline + intervention)
# ===================================================================
echo ""
echo "[6a] Copying multi-omics graph caches..."
CACHES_DST="$ZEN/graph_caches"
mkdir -p "$CACHES_DST"
# Human
cp "$ROOT/data/annotations/gtex_chr22_enhanced_cache.json" "$CACHES_DST/" 2>/dev/null \
    && echo "  ✓ human chr22 baseline cache" || echo "  ⚠ missing chr22 cache"
for f in "$ROOT/data/annotations/human_1kg_cache_chr"*.json; do
    [ -f "$f" ] && cp "$f" "$CACHES_DST/" && echo "  ✓ $(basename "$f")"
done
# Rice (12 chromosomes; v2 has heterogeneity prior_score)
for f in "$ROOT/data/rice_3k/annotations/rice_graph_cache_v2_Chr"*.json; do
    [ -f "$f" ] && cp "$f" "$CACHES_DST/" && echo "  ✓ $(basename "$f")"
done
for f in "$ROOT/data/rice_3k/annotations/rice_graph_cache_Chr"*.json; do
    [ -f "$f" ] && cp "$f" "$CACHES_DST/" && echo "  ✓ $(basename "$f") (baseline)"
done
# Arabidopsis (5 chromosomes; v3 has prior_score)
for f in "$ROOT/data/arabidopsis/arabidopsis_graph_cache_v3_chr"*.json \
         "$ROOT/data/arabidopsis/arabidopsis_graph_cache_chr"*.json; do
    [ -f "$f" ] && cp "$f" "$CACHES_DST/" && echo "  ✓ $(basename "$f")"
done
# Yeast
[ -f "$ROOT/data/yeast/yeast_graph_cache_v2.json" ] \
    && cp "$ROOT/data/yeast/yeast_graph_cache_v2.json" "$CACHES_DST/" \
    && echo "  ✓ yeast cache (intervention)"
[ -f "$ROOT/data/yeast/yeast_graph_cache.json" ] \
    && cp "$ROOT/data/yeast/yeast_graph_cache.json" "$CACHES_DST/" \
    && echo "  ✓ yeast cache (baseline)"

# ===================================================================
# 6b. IRRI 3kRG GWAS sumstats (the substrate for Sup Table S5)
# ===================================================================
echo ""
echo "[6b] Copying IRRI 3kRG GWAS sumstats..."
IRRI_SRC="$ROOT/data/rice_3k/results/irri_gwas"
IRRI_DST="$ZEN/irri_3krg_gwas"
mkdir -p "$IRRI_DST"
# Per-trait .glm.linear sumstats — gzip on the way in
for f in "$IRRI_SRC"/gwas_irri.*.glm.linear; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    if [ ! -f "$IRRI_DST/${base}.gz" ]; then
        gzip -c "$f" > "$IRRI_DST/${base}.gz" && echo "  ✓ ${base}.gz"
    fi
done
# Summary tables (already small) — copy if not staged
for f in irri_lambda_gc.tsv irri_gwas_summary.tsv irri_lead_loci.tsv; do
    [ -f "$IRRI_SRC/$f" ] && [ ! -f "$IRRI_DST/$f" ] \
        && cp "$IRRI_SRC/$f" "$IRRI_DST/" && echo "  ✓ $f"
done
# Phenotype summary + IRGC <-> 3kRG crosswalk
[ -f "$ROOT/data/rice_3k/pheno/irri_phenotype_summary.tsv" ] && [ ! -f "$IRRI_DST/irri_phenotype_summary.tsv" ] \
    && cp "$ROOT/data/rice_3k/pheno/irri_phenotype_summary.tsv" "$IRRI_DST/" \
    && echo "  ✓ irri_phenotype_summary.tsv"
[ -f "$ROOT/data/rice_3k/pheno/irri_irgc_to_3krg.tsv" ] \
    && cp "$ROOT/data/rice_3k/pheno/irri_irgc_to_3krg.tsv" "$IRRI_DST/" \
    && echo "  ✓ irri_irgc_to_3krg.tsv (sample-ID crosswalk)" \
    || echo "  ⚠ irri_irgc_to_3krg.tsv missing — please regenerate"
# Fine-mapping summary
for f in "$ROOT/data/rice_3k/results/irri_finemap_v2"/*summary*.tsv; do
    [ -f "$f" ] && cp "$f" "$IRRI_DST/" && echo "  ✓ $(basename "$f")"
done

# ===================================================================
# 6c. Cross-species lead lists (small TSVs already staged manually)
# ===================================================================
echo ""
echo "[6c] Refreshing cross-species lead lists..."
LEADS_SRC="$ROOT/results/multispecies_summary"
LEADS_DST="$ZEN/multispecies_leads"
mkdir -p "$LEADS_DST"
[ -f "$LEADS_SRC/cross_species_table.tsv" ] \
    && cp "$LEADS_SRC/cross_species_table.tsv" "$LEADS_DST/cross_species_summary.tsv" \
    && echo "  ✓ cross_species_summary.tsv"
[ -f "$LEADS_SRC/cross_species_report.md" ] \
    && cp "$LEADS_SRC/cross_species_report.md" "$LEADS_DST/" \
    && echo "  ✓ cross_species_report.md"
[ -f "$LEADS_SRC/heterogeneity_intervention_report.md" ] \
    && cp "$LEADS_SRC/heterogeneity_intervention_report.md" "$LEADS_DST/" \
    && echo "  ✓ heterogeneity_intervention_report.md"
for sp in arabidopsis rice yeast human human_gw; do
    src="$LEADS_SRC/validation_${sp}.tsv"
    case "$sp" in
        arabidopsis) dst="leads_arabidopsis_1001g.tsv" ;;
        rice)        dst="leads_rice_3krg.tsv" ;;
        yeast)       dst="leads_yeast_1011.tsv" ;;
        human)       dst="leads_human_chr22.tsv" ;;
        human_gw)    dst="leads_human_genome_wide.tsv" ;;
    esac
    [ -f "$src" ] && cp "$src" "$LEADS_DST/$dst" && echo "  ✓ $dst"
done

# ===================================================================
# 6d. Manuscript pre-print PDF
# ===================================================================
echo ""
echo "[6d] Copying manuscript PDF..."
[ -f "$ROOT/paper/finemapping_v1/main.pdf" ] \
    && cp "$ROOT/paper/finemapping_v1/main.pdf" "$ZEN/manuscript/manuscript_preprint.pdf" \
    && echo "  ✓ manuscript_preprint.pdf" \
    || echo "  ⚠ paper/finemapping_v1/main.pdf not built yet"

# ===================================================================
# 6e. LD references (genome-wide and chr22 ancestry-matched)
# ===================================================================
echo ""
echo "[6e] Copying LD reference matrices..."
LD_DST="$ZEN/ld_references"
mkdir -p "$LD_DST"
for f in "$ROOT/data/ld"/*ld.npz "$ROOT/data/1kg_NYGC_GRCh38"/*ld.npz; do
    [ -f "$f" ] && cp "$f" "$LD_DST/" && echo "  ✓ $(basename "$f")"
done
[ -f "$ROOT/data/panukb/panukb_ancestry_to_1kg_subset.tsv" ] \
    && cp "$ROOT/data/panukb/panukb_ancestry_to_1kg_subset.tsv" "$LD_DST/" \
    && echo "  ✓ panukb_ancestry_to_1kg_subset.tsv"

# ===================================================================
# 7. Source-code snapshot
# ===================================================================
echo ""
echo "[6/6] Creating source-code snapshot..."
cd "$ROOT"
SHA=$(git rev-parse HEAD 2>/dev/null || echo "unknown")
echo "$SHA" > "$ZEN/code_snapshot/commit_sha.txt"
# Prefer the latest semantic-version tag; fall back to HEAD
LATEST_TAG=$(git tag -l 'v[0-9]*' --sort=-v:refname | head -1)
if [ -n "$LATEST_TAG" ]; then
    SRC_REF="$LATEST_TAG"
    SRC_PREFIX="graphgwas-${LATEST_TAG#v}/"
    SRC_TGZ="graphgwas_${LATEST_TAG}.tar.gz"
else
    SRC_REF="HEAD"
    SRC_PREFIX="graphgwas-${SHA:0:8}/"
    SRC_TGZ="graphgwas_${SHA:0:8}.tar.gz"
fi
# Remove stale code-snapshot tarballs from older releases
find "$ZEN/code_snapshot" -name "graphgwas_v*.tar.gz" ! -name "$SRC_TGZ" -delete 2>/dev/null || true
git archive --format=tar.gz --prefix="$SRC_PREFIX" \
    -o "$ZEN/code_snapshot/$SRC_TGZ" "$SRC_REF" 2>/dev/null \
    && echo "  ✓ $SRC_TGZ (from $SRC_REF)" \
    || echo "  ⚠ git archive failed (not a git repo at this path)"

# ===================================================================
# Manifest with SHA-256
# ===================================================================
echo ""
echo "Computing SHA-256 checksums..."
{
    echo "# MANIFEST — GraphGWAS Zenodo deposit ${LATEST_TAG:-HEAD}"
    echo ""
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Project commit: $SHA"
    echo ""
    echo "| File | Size (bytes) | SHA-256 |"
    echo "|---|---:|---|"
    cd "$ZEN"
    find . -type f ! -name "MANIFEST.md" ! -name ".zenodo.json" -print0 \
    | sort -z \
    | while IFS= read -r -d '' f; do
        sz=$(stat -c %s "$f")
        sha=$(sha256sum "$f" | awk '{print $1}')
        echo "| \`$f\` | $sz | \`$sha\` |"
    done
} > "$ZEN/MANIFEST.md"
echo "  ✓ MANIFEST.md"

# ===================================================================
# 8. (Optional) Bundle Neo4j Community 5.26 tarball
#
# Set NEO4J_DOWNLOAD=1 to fetch ~160 MB from neo4j.com.
# Without this flag, the deposit ships dumps only and assumes users
# install Neo4j separately.
# ===================================================================
echo ""
echo "[8] Neo4j Community tarball (160 MB)..."
NEO4J_TGZ="$ZEN/graph_dumps/neo4j-community-5.26.0-unix.tar.gz"
if [ -f "$NEO4J_TGZ" ]; then
    echo "  ✓ already present at $NEO4J_TGZ"
elif [[ "${NEO4J_DOWNLOAD:-0}" == "1" ]]; then
    echo "  Downloading from dist.neo4j.org..."
    curl -L -o "$NEO4J_TGZ" \
        https://dist.neo4j.org/neo4j-community-5.26.0-unix.tar.gz \
        && echo "  ✓ neo4j-community-5.26.0-unix.tar.gz" \
        || echo "  ⚠ download failed; bundle without engine, document URL in README"
else
    echo "  ⏭ skipped (set NEO4J_DOWNLOAD=1 to fetch ~160 MB from dist.neo4j.org)"
fi

# ===================================================================
# 9. Flat upload bundle (GraphMana-style drag-and-drop into Zenodo)
#
# Produces zenodo/upload/ with one .tar.gz per top-level subfolder
# plus standalone PDFs and dumps. This is what you drag into the
# Zenodo "New upload" form.
# ===================================================================
echo ""
echo "[9] Building flat upload bundle..."
UP="$ZEN/upload"
rm -rf "$UP"
mkdir -p "$UP"

# Per-subfolder tarballs
for sub in benchmark_outputs panukb_results figure_source_data \
           simulation_seeds irri_3krg_gwas multispecies_leads \
           validation_catalogues ld_references graph_caches; do
    if [ -d "$ZEN/$sub" ] && [ -n "$(ls -A "$ZEN/$sub" 2>/dev/null)" ]; then
        tar -czf "$UP/graphgwas_${sub}.tar.gz" -C "$ZEN" "$sub"
        echo "  ✓ graphgwas_${sub}.tar.gz"
    fi
done

# Standalone large items (NOT re-tarballed — Zenodo serves them directly)
for f in graph_dumps/yeast_1011_v0.1.dump \
         graph_dumps/human_1kg_multiomics_v0.1.dump \
         graph_dumps/neo4j-community-5.26.0-unix.tar.gz \
         manuscript/manuscript_preprint.pdf \
         code_snapshot/$SRC_TGZ \
         CITATION.cff MANIFEST.md README.md; do
    if [ -f "$ZEN/$f" ]; then
        cp "$ZEN/$f" "$UP/$(basename "$f")"
        echo "  ✓ $(basename "$f")"
    fi
done

echo ""
echo "=== Staging complete ==="
echo "Subfolder tree (full hierarchy): $ZEN/"
echo "Flat upload bundle (drag into Zenodo): $UP/"
echo ""
echo "Bundle contents:"
ls -lh "$UP" | awk 'NR>1 {printf "  %-50s %s\n", $9, $5}'
echo ""
echo "Total upload bundle size:"
du -sh "$UP"
echo ""
echo "Next: open https://zenodo.org/uploads/new , drag every file in"
echo "      $UP into the form, then paste fields from $ZEN/.zenodo.json"
