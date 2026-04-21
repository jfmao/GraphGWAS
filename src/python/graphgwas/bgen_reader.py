"""BGEN-based genotype reader for hybrid GraphGWAS architecture.

For UKB-scale workloads, individual-level genotypes live in BGEN files rather
than Neo4j (storing 500K x 93M genotypes as edges is impractical). This module
provides a per-locus streaming interface matching the shape of Neo4j-backed
genotype queries: a variant table + a dense dosage matrix.

Usage:
    reader = BgenReader("/path/to/1kGP_bgen")
    variants, dosages = reader.load_locus("chr22", 16_050_000, 16_100_000)
    # variants: pd.DataFrame with columns [chr, pos, a1, a2]
    # dosages:  np.ndarray shape (n_samples, n_variants), float32, values in [0, 2]

Design notes:
* One BGEN per chromosome, named chr{1..22,X,Y}.bgen (PLINK2 output convention).
* Uses jeremymcrae/bgen under the hood. Needs no .bgi index: iterates by
  pre-loaded position array (binary search) and fetches variants by integer
  index for constant-time access.
* Dosage = P(het) + 2 * P(hom_alt), range [0, 2].
* Position array is cached per chromosome after first access (~1s for 1M variants).
"""
from __future__ import annotations

import bisect
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from bgen import BgenReader as _BgenFile


class BgenReader:
    """Per-chromosome BGEN reader. Caches file handles and position arrays."""

    def __init__(self, bgen_dir: Path | str, name_template: str = "chr{chr}.bgen"):
        self.dir = Path(bgen_dir)
        self.name_template = name_template
        self._files: dict[str, _BgenFile] = {}
        self._positions: dict[str, list[int]] = {}

    @staticmethod
    def _norm_chr(chrom: str) -> str:
        return chrom[3:] if chrom.startswith("chr") else chrom

    def _open(self, chrom: str) -> _BgenFile:
        stem = self._norm_chr(chrom)
        if stem not in self._files:
            path = self.dir / self.name_template.format(chr=stem)
            if not path.exists():
                raise FileNotFoundError(f"No BGEN for chromosome {chrom}: {path}")
            self._files[stem] = _BgenFile(str(path), delay_parsing=True)
        return self._files[stem]

    def _pos_array(self, chrom: str) -> list[int]:
        stem = self._norm_chr(chrom)
        if stem not in self._positions:
            bfile = self._open(chrom)
            self._positions[stem] = bfile.positions()
        return self._positions[stem]

    def n_samples(self, chrom: str = "22") -> int:
        return int(len(self._open(chrom).samples))

    def samples(self, chrom: str = "22") -> np.ndarray:
        return np.asarray(self._open(chrom).samples)

    def load_locus(
        self,
        chrom: str,
        start: int,
        end: int,
        format: Literal["dosage", "hardcall"] = "dosage",
    ) -> tuple[pd.DataFrame, np.ndarray]:
        """Return (variant_df, genotype_matrix) for a locus.

        variant_df columns: chr, pos, a1, a2
        genotype_matrix: shape (n_samples, n_variants)
          - format='dosage':   float32, values in [0, 2]
          - format='hardcall': int8, values in {0, 1, 2}
        """
        bfile = self._open(chrom)
        positions = self._pos_array(chrom)
        lo = bisect.bisect_left(positions, start)
        hi = bisect.bisect_right(positions, end)
        n = hi - lo
        if n == 0:
            empty = pd.DataFrame(columns=["chr", "pos", "a1", "a2"])
            return empty, np.empty((self.n_samples(chrom), 0), dtype=np.float32)
        n_samples = self.n_samples(chrom)
        dosage = np.empty((n_samples, n), dtype=np.float32)
        rows: list[tuple[str, int, str, str]] = []
        chr_label = f"chr{self._norm_chr(chrom)}"
        for j, idx in enumerate(range(lo, hi)):
            v = bfile[idx]
            # Dosage of the FIRST allele (allele index 0). plink2's VCF→BGEN conversion
            # places ALT as allele[0] and REF as allele[1], so dosage[0] = dosage of ALT.
            # Unphased biallelic (N,3): probs = (P(0/0), P(0/1), P(1/1))
            #                           dosage of a0 = 2*P(0/0) + P(0/1)
            # Phased biallelic (N,4):   probs = (P_mat(a0), P_mat(a1), P_pat(a0), P_pat(a1))
            #                           dosage of a0 = P_mat(a0) + P_pat(a0) = p[:,0] + p[:,2]
            p = v.probabilities
            if p.shape[1] == 4:
                dosage[:, j] = p[:, 0] + p[:, 2]
            elif p.shape[1] == 3:
                dosage[:, j] = 2.0 * p[:, 0] + p[:, 1]
            else:
                dosage[:, j] = v.minor_allele_dosage  # fallback
            a1, a2 = v.alleles[0], v.alleles[1] if len(v.alleles) > 1 else ""
            rows.append((chr_label, int(v.pos), a1, a2))
        variants = pd.DataFrame(rows, columns=["chr", "pos", "a1", "a2"])
        if format == "hardcall":
            return variants, np.rint(dosage).astype(np.int8)
        return variants, dosage

    def close(self) -> None:
        for f in self._files.values():
            try:
                f.close()
            except Exception:  # noqa: BLE001
                pass
        self._files.clear()
        self._positions.clear()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def _cli() -> None:
    import argparse, time
    ap = argparse.ArgumentParser(description="BGEN reader sanity checks")
    ap.add_argument("--bgen-dir", required=True, type=Path)
    ap.add_argument("--chr", default="22")
    ap.add_argument("--start", type=int, default=16_050_000)
    ap.add_argument("--end", type=int, default=16_100_000)
    args = ap.parse_args()
    reader = BgenReader(args.bgen_dir)
    t0 = time.time()
    variants, dosage = reader.load_locus(args.chr, args.start, args.end)
    dt = time.time() - t0
    print(f"Loaded {len(variants)} variants in {dt:.2f}s, dosage shape {dosage.shape}")
    print(variants.head())
    if len(variants):
        af = dosage.mean(axis=0) / 2.0
        print(f"AF range: {af.min():.4f}-{af.max():.4f}, mean {af.mean():.4f}")
        mac = dosage.sum(axis=0)
        print(f"MAC range: {mac.min():.0f}-{mac.max():.0f}")


if __name__ == "__main__":
    _cli()
