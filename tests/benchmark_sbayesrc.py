"""SBayesRC integration test (Wu et al. 2024 GWFM).

SBayesRC is a genome-wide fine-mapping method (multicomponent Bayesian mixture
with functional annotations) that operates on summary statistics for all
HapMap3 SNPs simultaneously. It's the central new method in Zheng et al.
2024 Nature Genetics and Wu et al. 2026 Nature Genetics (Fig. 4b benchmark).

Architecture mismatch: SBayesRC fits ~1.2M HM3 SNPs in one MCMC pass,
optimised for biobank-scale (N ≥ 100K). Our F1 single-causal-in-LD
simulations on 1KG (N=3,202) are designed for region-specific Bayesian
methods (SuSiE, FINEMAP, HBP). A direct head-to-head on our F1 setup
under-uses SBayesRC.

This script does two things:
1. Verifies the SBayesRC pipeline runs end-to-end on the canonical example
   (Wu lab UKB sumstats, N≈340K, full HM3 LD ref).
2. Documents output format and timing for the paper.

Output: results/benchmark_v2/sbayesrc_integration/sbayesrc_summary.{json,md}
"""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

OUT = Path("/mnt/data/GraphGWAS/results/benchmark_v2/sbayesrc_integration")
OUT.mkdir(parents=True, exist_ok=True)
SBAYES_DIR = Path("/mnt/data/GraphGWAS/data/sbayesrc")


def run_pipeline():
    """Run tidy → sbayesrc on the example UKB sumstats."""
    log: dict[str, dict] = {}

    # Step 1: tidy
    tidy_script = """
    library(SBayesRC)
    SBayesRC::tidy(mafile="example.ma", LDdir="ukbEUR_HM3",
                   output="example_tidy.ma", log2file=FALSE)
    """
    t0 = time.time()
    proc = subprocess.run(
        ["R", "--no-save", "-e", tidy_script],
        cwd=str(SBAYES_DIR), capture_output=True, text=True, timeout=600,
    )
    log["tidy"] = {
        "runtime_s": time.time() - t0,
        "exit_code": proc.returncode,
        "stdout_tail": proc.stdout.strip().split("\n")[-3:],
    }

    # Step 2: sbayesrc with annotations
    main_script = """
    library(SBayesRC)
    SBayesRC::sbayesrc(mafile="example_tidy.ma", LDdir="ukbEUR_HM3",
                       outPrefix="example_sbrc",
                       annot="annot_baseline2.2.txt",
                       niter=500, burn=200, log2file=FALSE)
    """
    t0 = time.time()
    proc = subprocess.run(
        ["R", "--no-save", "-e", main_script],
        cwd=str(SBAYES_DIR), capture_output=True, text=True, timeout=3600,
    )
    log["sbayesrc"] = {
        "runtime_s": time.time() - t0,
        "exit_code": proc.returncode,
        "stdout_tail": proc.stdout.strip().split("\n")[-15:],
    }

    # Step 3: parse output
    out_file = SBAYES_DIR / "example_sbrc.txt"
    if out_file.exists():
        # Parse: SNP A1 BETA SE PIP BETAlast
        lines = out_file.read_text().split("\n")
        n_snps = len(lines) - 2  # header + trailing newline
        # Top PIPs
        rows = []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    rows.append((parts[0], float(parts[2]), float(parts[4])))
                except ValueError:
                    pass
        rows.sort(key=lambda r: -r[2])
        top_pip_snps = rows[:20]
        log["output"] = {
            "n_snps": n_snps,
            "n_pip_gt_0_5": sum(1 for r in rows if r[2] > 0.5),
            "n_pip_gt_0_9": sum(1 for r in rows if r[2] > 0.9),
            "max_pip": rows[0][2] if rows else None,
            "top_20": [{"snp": r[0], "beta": r[1], "pip": r[2]}
                       for r in top_pip_snps],
        }
    return log


def main():
    print("Running SBayesRC pipeline on example UKB sumstats...")
    log = run_pipeline()
    print(f"\n  tidy:     {log['tidy']['runtime_s']:.1f}s, exit={log['tidy']['exit_code']}")
    print(f"  sbayesrc: {log['sbayesrc']['runtime_s']:.1f}s, exit={log['sbayesrc']['exit_code']}")
    if "output" in log:
        o = log["output"]
        print(f"  output:   {o['n_snps']:,} SNPs, "
              f"{o['n_pip_gt_0_5']:,} with PIP>0.5, "
              f"{o['n_pip_gt_0_9']:,} with PIP>0.9, max PIP {o['max_pip']:.3f}")
        print("\nTop 5 PIPs:")
        for r in o["top_20"][:5]:
            print(f"  {r['snp']:20s} β={r['beta']:+.3f} PIP={r['pip']:.4f}")

    (OUT / "sbayesrc_summary.json").write_text(json.dumps(log, indent=2, default=str))

    md = []
    md.append("# SBayesRC Integration Test\n")
    md.append(f"- tidy runtime: **{log['tidy']['runtime_s']:.1f} s**")
    md.append(f"- sbayesrc runtime (500 MCMC iter): **{log['sbayesrc']['runtime_s']:.1f} s**")
    if "output" in log:
        md.append(f"- output SNPs: **{log['output']['n_snps']:,}**")
        md.append(f"- SNPs with PIP > 0.5: **{log['output']['n_pip_gt_0_5']:,}**")
        md.append(f"- SNPs with PIP > 0.9: **{log['output']['n_pip_gt_0_9']:,}**")
        md.append(f"- max PIP: **{log['output']['max_pip']:.4f}**\n")
        md.append("## Top 20 PIPs:\n")
        md.append("| SNP | β | PIP |")
        md.append("|---|---:|---:|")
        for r in log["output"]["top_20"]:
            md.append(f"| {r['snp']} | {r['beta']:+.4f} | {r['pip']:.4f} |")
    (OUT / "sbayesrc_summary.md").write_text("\n".join(md))
    print(f"\nWrote {OUT}/sbayesrc_summary.{{json,md}}")


if __name__ == "__main__":
    main()
