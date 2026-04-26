"""Pan-UKB data client for summary-statistics-only fine-mapping.

Fetches per-locus summary statistics and in-sample LD matrices from the
public Pan-UK Biobank release on Amazon S3:

    s3://pan-ukb-us-east-1/sumstats_release/*.tsv.bgz
    s3://pan-ukb-us-east-1/ld_release/UKBB.{ANC}.*

No AWS credentials required — the bucket is fully public. Hail is an
optional dependency; only the LD-slicing functions need it. The sumstats
fetcher uses pandas over HTTPS and has no Hail dependency.

Coordinate system: Pan-UKB sumstats are GRCh37. The LD variant manifest
ships in both GRCh37 (`UKBB.{ANC}.ldadj.variant.ht`) and GRCh38
(`UKBB.{ANC}.ldadj.variant.b38.ht`). We default to GRCh38 so annotations
join directly against our GENCODE v47 / GTEx v8 graph.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
import io
import shutil
import subprocess
import urllib.request
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd


# ===================================================================
# Constants
# ===================================================================

PANUKB_S3_HTTPS = "https://pan-ukb-us-east-1.s3.amazonaws.com"
PANUKB_S3 = "s3://pan-ukb-us-east-1"

ANCESTRIES = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")

# Sample sizes from the Pan-UKB manuscript (Karczewski et al. 2024).
# These are approximate per-ancestry cohort sizes used for LD matrix
# computation and as the default meta-analysis N.
ANCESTRY_N: dict[str, int] = {
    "EUR": 420_531,
    "CSA": 8_876,
    "AFR": 6_636,
    "EAS": 2_709,
    "MID": 1_599,
    "AMR": 980,
}

# Ancestries we recommend using for fine-mapping (N > 2500). AMR and
# MID are available but too small for reliable single-ancestry LD
# estimation.
RECOMMENDED_ANCESTRIES = ("EUR", "CSA", "AFR", "EAS")

PHENO_MANIFEST_URL = (
    f"{PANUKB_S3_HTTPS}/sumstats_release/phenotype_manifest.tsv.bgz"
)
VARIANT_MANIFEST_URL = (
    f"{PANUKB_S3_HTTPS}/sumstats_release/full_variant_qc_metrics.txt.bgz"
)


# ===================================================================
# Data classes
# ===================================================================

@dataclass
class LocusSumstats:
    """Per-locus per-ancestry Pan-UKB summary statistics.

    Attributes:
        chr, start, end: GRCh37 locus window.
        ancestry: one of ANCESTRIES.
        variants: DataFrame with columns
            variant_id, chr, pos, ref, alt, af, beta, se, z, log10p,
            low_confidence.
        n_samples: effective N for this ancestry.
    """

    chr: str
    start: int
    end: int
    ancestry: str
    variants: pd.DataFrame
    n_samples: int


@dataclass
class LDSlice:
    """Per-locus LD matrix slice from Pan-UKB.

    Attributes:
        chr, start, end: locus window.
        ancestry: one of ANCESTRIES.
        variant_ids: list of variant IDs (chr:pos:ref:alt), GRCh38 by default.
        R: n × n signed-correlation matrix (float32). `R**2` gives r².
        build: "b38" or "b37".
    """

    chr: str
    start: int
    end: int
    ancestry: str
    variant_ids: list[str]
    R: np.ndarray
    build: str = "b38"


# ===================================================================
# Sumstats fetcher (no Hail dependency)
# ===================================================================

def _read_bgz_tsv(url: str, **read_csv_kwargs) -> pd.DataFrame:
    """Fetch a bgzipped TSV over HTTPS and return a pandas DataFrame.

    Pan-UKB ships all TSVs as .tsv.bgz (block-gzip). pandas handles
    standard gzip transparently; bgzip is a gzip superset and loads
    identically unless index-based random access is needed.
    """
    if not url.startswith(("https://", "http://")):
        raise ValueError(
            f"Refusing non-HTTP URL: {url!r}. "
            "Pan-UKB fetcher only accepts HTTP(S) endpoints."
        )
    with urllib.request.urlopen(url) as resp:  # nosec B310 — scheme guarded above
        buf = io.BytesIO(resp.read())
    return pd.read_csv(
        buf, sep="\t", compression="gzip", **read_csv_kwargs,
    )


def load_phenotype_manifest() -> pd.DataFrame:
    """Load the Pan-UKB phenotype manifest (all ~7,228 phenotypes).

    Columns include: trait_type, phenocode, pheno_sex, coding, modifier,
    description, n_cases_full_cohort_both_sexes, and per-ancestry
    sample counts, plus download URLs for the per-phenotype TSV files.
    """
    return _read_bgz_tsv(PHENO_MANIFEST_URL)


def _build_sumstats_fname(
    phenocode: str,
    trait_type: str,
    coding: str | None,
    modifier: str | None,
) -> str:
    """Compose Pan-UKB per-phenotype filename from manifest key parts.

    Pan-UKB filenames follow the pattern
    `{trait_type}-{phenocode}[-{coding}]-both_sexes[-{modifier}].tsv.bgz`
    with small per-trait-type variations. Continuous biomarkers are
    invariably inverse-rank-normal transformed and the file carries the
    `-irnt` modifier; categoricals add a coding. Always confirm via the
    manifest's `filename` column when possible.
    """
    parts = [trait_type, phenocode]
    if coding:
        parts.append(coding)
    parts.append("both_sexes")
    if modifier:
        parts.append(modifier)
    return "-".join(parts) + ".tsv.bgz"


_PANUKB_SUMSTATS_COLS = [
    "chr", "pos", "ref", "alt",
    "af_meta_hq", "beta_meta_hq", "se_meta_hq",
    "neglog10_pval_meta_hq", "neglog10_pval_heterogeneity_hq",
    "af_meta", "beta_meta", "se_meta",
    "neglog10_pval_meta", "neglog10_pval_heterogeneity",
    "af_AFR", "af_AMR", "af_CSA", "af_EAS", "af_EUR", "af_MID",
    "beta_AFR", "beta_AMR", "beta_CSA", "beta_EAS", "beta_EUR", "beta_MID",
    "se_AFR", "se_AMR", "se_CSA", "se_EAS", "se_EUR", "se_MID",
    "neglog10_pval_AFR", "neglog10_pval_AMR", "neglog10_pval_CSA",
    "neglog10_pval_EAS", "neglog10_pval_EUR", "neglog10_pval_MID",
    "low_confidence_AFR", "low_confidence_AMR", "low_confidence_CSA",
    "low_confidence_EAS", "low_confidence_EUR", "low_confidence_MID",
]


def fetch_sumstats_locus(
    phenocode: str,
    chr: str,
    start: int,
    end: int,
    ancestries: Iterable[str] = RECOMMENDED_ANCESTRIES,
    trait_type: str = "continuous",
    coding: str | None = None,
    modifier: str | None = None,
    tabix_bin: str | None = None,
) -> dict[str, LocusSumstats]:
    """Fetch per-ancestry sumstats for a locus from Pan-UKB via tabix.

    Uses `tabix` (htslib) over HTTPS to slice a ~100 kb–1 Mb locus out
    of the remote bgzipped TSV without downloading the full 2.3 GB
    file — typically ~1–5 s per query. Requires htslib's tabix binary
    on PATH (or supplied via `tabix_bin`).

    Args:
        phenocode: Pan-UKB phenocode (e.g. "21001" for BMI).
        chr: chromosome, GRCh37 ("1"..."22", "X").
        start, end: window (bp), GRCh37.
        ancestries: subset of ANCESTRIES.
        trait_type: "continuous", "biomarkers", "categorical",
            "icd10", "phecode", or "prescriptions".
        coding, modifier: additional manifest keys where needed
            (continuous/biomarkers use modifier="irnt").
        tabix_bin: path to tabix (defaults to $PATH lookup).

    Returns:
        dict mapping ancestry → LocusSumstats.
    """
    tabix = tabix_bin or shutil.which("tabix")
    if not tabix:
        raise RuntimeError(
            "tabix (htslib) not found on PATH. Install with "
            "`conda install -c bioconda htslib` or `apt install tabix`.",
        )

    fname = _build_sumstats_fname(phenocode, trait_type, coding, modifier)
    url = f"{PANUKB_S3_HTTPS}/sumstats_flat_files/{fname}"
    region = f"{chr}:{int(start)}-{int(end)}"
    res = subprocess.run(
        [tabix, url, region], capture_output=True, text=True, timeout=300,
    )
    if res.returncode != 0:
        raise RuntimeError(
            f"tabix query failed for {url} region {region}: {res.stderr.strip()}",
        )
    if not res.stdout:
        return {}

    df = pd.read_csv(
        io.StringIO(res.stdout), sep="\t", header=None,
        names=_PANUKB_SUMSTATS_COLS, dtype={"chr": str},
    )

    out: dict[str, LocusSumstats] = {}
    for anc in ancestries:
        beta_col = f"beta_{anc}"
        se_col = f"se_{anc}"
        if beta_col not in df.columns or se_col not in df.columns:
            continue
        sub = df[["chr", "pos", "ref", "alt",
                  beta_col, se_col, f"af_{anc}",
                  f"neglog10_pval_{anc}", f"low_confidence_{anc}"]].copy()
        sub.columns = ["chr", "pos", "ref", "alt",
                       "beta", "se", "af", "log10p", "low_confidence"]
        sub = sub.replace({"NA": pd.NA}).dropna(subset=["beta", "se"])
        if sub.empty:
            continue
        sub["beta"] = pd.to_numeric(sub["beta"], errors="coerce")
        sub["se"] = pd.to_numeric(sub["se"], errors="coerce")
        sub["af"] = pd.to_numeric(sub["af"], errors="coerce")
        sub["log10p"] = pd.to_numeric(sub["log10p"], errors="coerce")
        sub = sub.dropna(subset=["beta", "se"])
        if sub.empty:
            continue
        sub["variant_id"] = (sub["chr"].astype(str) + ":" + sub["pos"].astype(str)
                             + ":" + sub["ref"] + ":" + sub["alt"])
        sub["z"] = sub["beta"] / sub["se"]
        sub = sub.reset_index(drop=True)
        out[anc] = LocusSumstats(
            chr=str(chr), start=int(start), end=int(end), ancestry=anc,
            variants=sub, n_samples=ANCESTRY_N.get(anc, 0),
        )
    return out


# ===================================================================
# LD BlockMatrix slicing (Hail-dependent)
# ===================================================================

def _require_hail():
    """Import Hail on demand with a helpful error message."""
    try:
        import hail as hl  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "Pan-UKB LD slicing requires Hail >= 0.2.42. Install with:\n"
            "  pip install hail  (needs Java 11+, Spark)\n"
            "  or use a Hail-compatible conda environment.\n"
            "See https://hail.is/docs/0.2/getting_started.html"
        ) from e


def fetch_ld_slice(
    ancestry: str,
    chr: str,
    start: int,
    end: int,
    build: str = "b38",
    init_hail: bool = True,
) -> LDSlice:
    """Fetch an LD correlation slice for a locus from Pan-UKB.

    Loads the bias-adjusted Pan-UKB LD BlockMatrix for the requested
    ancestry, slices it to the variants in the window, and returns a
    dense correlation matrix. This function requires Hail.

    Storage note: the full BlockMatrix is 2.8–15.5 TB per ancestry,
    but Hail reads only the blocks needed for the slice (~1–3 GB per
    locus). Ensure your environment has read access to
    s3://pan-ukb-us-east-1/ (public, no auth).

    Args:
        ancestry: one of ANCESTRIES.
        chr: chromosome ("1"..."22", "X").
        start, end: window (bp) in the requested build.
        build: "b38" (default, joins our GENCODE v47 graph directly)
            or "b37" (matches Pan-UKB sumstats).
        init_hail: call hl.init() if no Hail context is active.

    Returns:
        LDSlice with signed-correlation R matrix and variant IDs.
    """
    assert ancestry in ANCESTRIES, f"Unknown ancestry: {ancestry}"
    assert build in ("b37", "b38"), "build must be b37 or b38"

    _require_hail()
    import hail as hl
    from hail.linalg import BlockMatrix

    if init_hail and not hl.utils.java.Env._hc:
        import tempfile
        hail_log = str(Path(tempfile.gettempdir()) / "hail.log")
        hl.init(default_reference="GRCh37" if build == "b37" else "GRCh38",
                log=hail_log, quiet=True)

    variant_suffix = "ldadj.variant.ht" if build == "b37" else "ldadj.variant.b38.ht"
    variant_path = f"{PANUKB_S3}/ld_release/UKBB.{ancestry}.{variant_suffix}"
    bm_path = f"{PANUKB_S3}/ld_release/UKBB.{ancestry}.ldadj.bm"

    # Read variant table, filter to window
    ht = hl.read_table(variant_path)
    ref_str = "GRCh37" if build == "b37" else "GRCh38"
    chr_str = str(chr) if build == "b37" else f"chr{chr}"
    interval = hl.parse_locus_interval(
        f"{chr_str}:{start}-{end}", reference_genome=ref_str,
    )
    ht = ht.filter(interval.contains(ht.locus))
    ht = ht.add_index("_idx")
    ht_local = ht.collect()

    if not ht_local:
        return LDSlice(chr=str(chr), start=start, end=end, ancestry=ancestry,
                       variant_ids=[], R=np.zeros((0, 0)), build=build)

    bm = BlockMatrix.read(bm_path)
    indices = [int(row["idx"]) for row in ht_local]
    sub = bm.filter(indices, indices).to_numpy()  # dense n × n

    variant_ids = [
        f"{row.locus.contig.lstrip('chr')}:{row.locus.position}:{row.alleles[0]}:{row.alleles[1]}"
        for row in ht_local
    ]
    return LDSlice(
        chr=str(chr), start=start, end=end, ancestry=ancestry,
        variant_ids=variant_ids, R=sub.astype(np.float32), build=build,
    )


# ===================================================================
# Sumstats + LD alignment
# ===================================================================

def align_sumstats_to_ld(
    sumstats: LocusSumstats,
    ld: LDSlice,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Intersect sumstats with LD variants, aligning order and signs.

    Returns (z, R, variant_ids) for the intersection, with z signs
    flipped where sumstats and LD use opposite reference alleles.
    """
    ss = sumstats.variants.set_index("variant_id")
    # Build LD variant index → position for fast lookup
    ld_idx_by_vid = {vid: i for i, vid in enumerate(ld.variant_ids)}

    keep_ld_idx: list[int] = []
    keep_ss_idx: list[int] = []
    flips: list[int] = []
    for pos, vid in enumerate(ld.variant_ids):
        if vid in ss.index:
            keep_ld_idx.append(pos)
            keep_ss_idx.append(ss.index.get_loc(vid))
            flips.append(1)
            continue
        # Try allele swap: chr:pos:alt:ref
        chr_p, p, a, b = vid.split(":")
        alt_vid = f"{chr_p}:{p}:{b}:{a}"
        if alt_vid in ss.index:
            keep_ld_idx.append(pos)
            keep_ss_idx.append(ss.index.get_loc(alt_vid))
            flips.append(-1)

    if not keep_ld_idx:
        return np.zeros(0), np.zeros((0, 0)), []

    ld_idx = np.asarray(keep_ld_idx, dtype=int)
    ss_vals = sumstats.variants.iloc[keep_ss_idx]
    flip = np.asarray(flips, dtype=np.float32)
    z = ss_vals["z"].to_numpy(dtype=np.float32) * flip
    R = ld.R[np.ix_(ld_idx, ld_idx)]
    # Apply sign flips symmetrically to R (if variant i is flipped, its
    # row and column correlations must flip too)
    R = R * flip[:, None] * flip[None, :]
    aligned_vids = [ld.variant_ids[i] for i in keep_ld_idx]
    return z, R, aligned_vids
