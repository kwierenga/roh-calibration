"""
24_bootstrap_lstar.py - bootstrap confidence intervals on the decisive ROH
length L* from the trio-children background null (script 21).

WHY: script 21 reports a POINT estimate L*_screened ~1.55-1.65 Mb per
superpopulation. Reviewers will ask how stable that number is. This resamples
the natural sampling unit -- the trio CHILD (a nominally-outbred child's two
homologs are one draw of the population's chance/background autozygosity) --
with replacement, recomputes L*_screened on each bootstrap replicate, and
reports a percentile confidence interval.

Reuses script 21's single VCF parse (parse_children) and its empirical-null
aggregator (agg_null), so the point estimate here equals 21's by construction;
only the resampling is added. One parse, B cheap recomputes. Children flagged
as F_ROH outliers (recent shared ancestry) are excluded before resampling, the
same screen 21 applies to L*_screened.

Outputs:
  trio_null_lstar_bootstrap.tsv   per-pop point, CI, median, SE, non-converge
  trio_null_lstar_bootstrap.txt   human-readable summary

Usage:
  python 24_bootstrap_lstar.py chr22            # smoke test (fast)
  python 24_bootstrap_lstar.py                  # all 22 autosomes (background)
  python 24_bootstrap_lstar.py --B=2000 chr22   # override replicate count
"""

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent

# script 21's module name starts with a digit -> load via importlib
_spec = importlib.util.spec_from_file_location(
    "trionull21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m21)

POPS = m21.POPS
OUTLIER_F = m21.OUTLIER_F
L_GRID = m21.L_GRID

B_DEFAULT = 1000
RNG_SEED = 17
CI = (2.5, 97.5)
# replicates that never cross the decision threshold get this ">grid" sentinel
INF_SENTINEL = float(L_GRID[-1]) + 0.05

OUT_TSV = HERE / "trio_null_lstar_bootstrap.tsv"
OUT_TXT = HERE / "trio_null_lstar_bootstrap.txt"


def main():
    t0 = time.time()
    args = sys.argv[1:]
    bflag = [a for a in args if a.startswith("--B=")]
    B = int(bflag[0].split("=", 1)[1]) if bflag else B_DEFAULT
    args = [a for a in args if not a.startswith("--B=")]
    chroms = args or [f"chr{n}" for n in range(1, 23)]

    print(f"parsing children (chroms={chroms}) ...", flush=True)
    childsegs, childtot, childspan, nchild = m21.parse_children(chroms)

    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for p in POPS:
        if childsegs[p] is None:
            continue
        n = nchild[p]
        froh = np.array([childtot[p][j] / childspan[p][j] if childspan[p][j] else 0.0
                         for j in range(n)])
        keep = np.flatnonzero(froh <= OUTLIER_F)
        m = keep.size
        if m == 0:
            continue

        # point estimate -- identical to script 21's L*_screened by construction
        _, Lpoint = m21.agg_null([childsegs[p][j] for j in keep],
                                 [childspan[p][j] for j in keep])

        boot = np.empty(B)
        n_inf = 0
        for b in range(B):
            idx = keep[rng.integers(0, m, size=m)]   # resample children w/ replacement
            _, Lb = m21.agg_null([childsegs[p][j] for j in idx],
                                 [childspan[p][j] for j in idx])
            if not np.isfinite(Lb):
                n_inf += 1
                Lb = INF_SENTINEL
            boot[b] = Lb

        lo, hi = np.percentile(boot, CI)
        rows.append({
            "pop": p, "n_screened": int(m), "n_outliers": int(n - m),
            "Lpoint": float(Lpoint), "ci_lo": float(lo), "ci_hi": float(hi),
            "median": float(np.median(boot)), "se": float(np.std(boot, ddof=1)),
            "n_nonconv": n_inf,
        })
        print(f"  {p}: L*={Lpoint:.2f}  95% CI [{lo:.2f}, {hi:.2f}]  "
              f"(n_screened={m}, nonconv={n_inf})", flush=True)

    with OUT_TSV.open("w", encoding="utf-8") as fh:
        fh.write("pop\tn_screened\tn_outliers\tLstar_point_Mb\tCI2.5_Mb\tCI97.5_Mb\t"
                 "boot_median_Mb\tboot_SE_Mb\tn_nonconverge\tB\n")
        for r in rows:
            fh.write(f"{r['pop']}\t{r['n_screened']}\t{r['n_outliers']}\t"
                     f"{r['Lpoint']:.2f}\t{r['ci_lo']:.2f}\t{r['ci_hi']:.2f}\t"
                     f"{r['median']:.2f}\t{r['se']:.3f}\t{r['n_nonconv']}\t{B}\n")

    lines = [
        "# Bootstrap CI on decisive ROH length L* (trio-children background null)",
        f"# chroms={','.join(chroms)}  B={B}  resampling unit = trio child (screened)",
        f"# OUTLIER_F={OUTLIER_F}  CI={CI[0]}-{CI[1]} percentile  seed={RNG_SEED}",
        f"# non-converging replicates set to L*>{L_GRID[-1]:.2f} Mb (sentinel {INF_SENTINEL:.2f})",
        f"# wall={time.time()-t0:.0f}s\n",
        "pop\tn_scr\tL*_point\t95% CI (Mb)\tboot_median\tSE\tn_nonconv",
    ]
    for r in rows:
        lines.append(f"{r['pop']}\t{r['n_screened']}\t{r['Lpoint']:.2f}\t"
                     f"[{r['ci_lo']:.2f}, {r['ci_hi']:.2f}]\t{r['median']:.2f}\t"
                     f"{r['se']:.3f}\t{r['n_nonconv']}")
    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}")


if __name__ == "__main__":
    main()
