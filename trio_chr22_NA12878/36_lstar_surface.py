"""
36_lstar_surface.py - the decisive ROH length L* as a SURFACE over the two
analyst choices it is most sensitive to, with bootstrap CIs.

WHY (reviewer hardening, item 2): script 21 reports a single point estimate
L*_screened ~1.55-1.65 Mb at GAP_TOL=1 and OUTLIER_F=0.0156. Two objections:
  (a) the headline depends on the gap tolerance -- at gap=0 (strict, no tolerated
      mismatch) L* collapses toward the analytic ~0.7-0.85 Mb; at gap=2 it rises
      toward ~2.1-2.6 Mb. The single number hides that dependence.
  (b) the cryptic-relatedness screen (F_ROH > OUTLIER_F) removes a large fraction
      of some populations (notably AMR), so "ancestry-robust convergence" is partly
      a product of aggressive screening.
This script makes both explicit instead of burying them:
  TABLE A  L* over gap in {0,1,2} x screen in {all, screened}, per population,
           each with a bootstrap 95% CI (resampling unit = trio child).
  TABLE B  at the operating gap=1, L* and the RETAINED sample size as the screen
           cutoff OUTLIER_F is swept {none, 0.05, 0.0156, 0.005}, so the reader
           sees exactly how much L* moves with the screen and how many children
           survive it (the AMR-screening criticism, quantified).

The expensive VCF I/O is done ONCE; segments for all three gap tolerances are
accumulated in that single pass. The screen (F_ROH) is computed at a FIXED gap
(SCREEN_GAP=1) and applied across the gap axis, so the gap axis isolates the
run-merging effect and the screen axis isolates the relatedness-screen effect.

Reuses script 21's VCF constants and its empirical-null aggregator agg_null
(gap-independent), so a screened gap=1 cell reproduces script 21's L*_screened.

Outputs:
  lstar_surface.tsv   machine-readable surface (Table A) + screen sweep (Table B)
  lstar_surface.txt   human-readable
Usage:
  python 36_lstar_surface.py chr22          # smoke test (fast)
  python 36_lstar_surface.py                # all 22 autosomes (background)
  python 36_lstar_surface.py --B=2000 chr22 # override bootstrap replicate count
"""

import gzip
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent

# script 21's module name starts with a digit -> load via importlib, reuse its
# constants + agg_null so a screened gap=1 cell matches 21 by construction.
_spec = importlib.util.spec_from_file_location(
    "trionull21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m21)

POPS = m21.POPS
MAF_MIN = m21.MAF_MIN
MAX_SNP_GAP_BP = m21.MAX_SNP_GAP_BP
MIN_KEEP_MB = m21.MIN_KEEP_MB
FROH_MIN_MB = m21.FROH_MIN_MB
AF_PRE = m21.AF_PRE
DATA_DIR_OTHER = m21.DATA_DIR_OTHER
agg_null = m21.agg_null          # uses THR_PC, L_GRID, C_IBD (all gap-independent)
L_GRID = m21.L_GRID

GAPS = [0, 1, 2]
SCREEN_GAP = 1                   # gap at which F_ROH is computed for the screen
OUTLIER_LADDER = [None, 0.05, 0.0156, 0.005]   # screen-cutoff sweep (Table B)
OPERATING_GAP = 1
B_DEFAULT = 1000
RNG_SEED = 17
CI = (2.5, 97.5)
INF_SENTINEL = float(L_GRID[-1]) + 0.05

OUT_TSV = HERE / "lstar_surface.tsv"
OUT_TXT = HERE / "lstar_surface.txt"


def roh_lengths_g(hom, pos, gap_tol):
    """ROH segment lengths (Mb) with a given run-merging gap tolerance. Identical
    logic to script 21's roh_lengths but with gap_tol as a parameter."""
    m = hom.copy()
    if gap_tol > 0:
        pad = np.concatenate(([1], m.astype(np.int8), [1])); dd = np.diff(pad)
        hs = np.flatnonzero(dd == -1); he = np.flatnonzero(dd == 1)
        short = (he - hs) <= gap_tol
        if short.any():
            diff = np.zeros(m.size + 1, dtype=np.int32)
            np.add.at(diff, hs[short], 1); np.add.at(diff, he[short], -1)
            m = m | (np.cumsum(diff[:-1]) > 0)
    n = m.size
    intra = np.zeros(n, bool)
    intra[1:] = m[1:] & m[:-1] & ((pos[1:] - pos[:-1]) <= MAX_SNP_GAP_BP)
    starts = np.flatnonzero(m & ~intra)
    ends = m.copy(); ends[:-1] &= ~intra[1:]; ends = np.flatnonzero(ends)
    if starts.size == 0:
        return np.empty(0, np.float32)
    return ((pos[ends] - pos[starts]) / 1e6).astype(np.float32)


def parse_multigap(chroms):
    """One VCF pass; accumulate per-child ROH segments for every gap in GAPS, plus
    F_ROH burden (at SCREEN_GAP) and scanned span for the screen.

    Returns segs_by_gap[gap][p][j] (list of per-chrom seg arrays), childtot[p][j]
    (ROH>=FROH_MIN_MB burden Mb at SCREEN_GAP), childspan[p][j] (Mb), nchild[p]."""
    t0 = time.time()
    sp = m21.load_superpop()
    kids = m21.load_children(sp)
    segs_by_gap = {g: {p: None for p in POPS} for g in GAPS}
    childtot = {p: None for p in POPS}
    childspan = {p: None for p in POPS}
    nchild = {p: 0 for p in POPS}

    for chrom in chroms:
        vcf = HERE / "chr22_phased.vcf.gz" if chrom == "chr22" else DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz"
        if not vcf.exists():
            print(f"  [{chrom}] SKIP"); continue
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]; break
        colpop = {p: [] for p in POPS}
        for j, s in enumerate(samples):
            if s in kids:
                colpop[kids[s]].append(j)
        for p in POPS:
            if childtot[p] is None and colpop[p]:
                n = len(colpop[p])
                for g in GAPS:
                    segs_by_gap[g][p] = [[] for _ in range(n)]
                childtot[p] = [0.0] * n; childspan[p] = [0.0] * n; nchild[p] = n
        rows = {p: [] for p in POPS}; pos = {p: [] for p in POPS}
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
                    continue
                common = []
                for kv in f[7].split(";"):
                    if kv[:3] != "AF_":
                        continue
                    for p, pre in AF_PRE.items():
                        if kv.startswith(pre):
                            try:
                                af = float(kv[len(pre):])
                            except ValueError:
                                af = None
                            if af is not None and min(af, 1 - af) >= MAF_MIN:
                                common.append(p)
                            break
                if not common:
                    continue
                gts = f[9:]; p1 = int(f[1])
                for p in common:
                    cols = colpop[p]
                    if not cols:
                        continue
                    rows[p].append(bytes(1 if gts[c][0] == gts[c][2] else 0 for c in cols))
                    pos[p].append(p1)
        for p in POPS:
            if not rows[p]:
                continue
            ncol = len(colpop[p])
            mat = np.frombuffer(b"".join(rows[p]), dtype=np.int8).reshape(len(rows[p]), ncol)
            pa = np.asarray(pos[p], dtype=np.int64)
            span = (pa[-1] - pa[0]) / 1e6
            for j in range(ncol):
                hom = mat[:, j].astype(bool)
                childspan[p][j] += span
                for g in GAPS:
                    sl = roh_lengths_g(hom, pa, g)
                    slk = sl[sl > MIN_KEEP_MB]
                    if slk.size:
                        segs_by_gap[g][p][j].append(slk)
                    if g == SCREEN_GAP:
                        childtot[p][j] += float(sl[sl >= FROH_MIN_MB].sum())
        print(f"  [{chrom}] children/pop " + " ".join(f"{p}:{nchild[p]}" for p in POPS)
              + f"  ({time.time()-t0:.0f}s)")
        sys.stdout.flush()
    return segs_by_gap, childtot, childspan, nchild


def lstar_ci(segs_p, span_p, keep_idx, rng, B):
    """Point L* over the kept children + a bootstrap percentile CI (resample
    children with replacement). Returns (point, lo, hi, n_nonconv)."""
    _, point = agg_null([segs_p[j] for j in keep_idx], [span_p[j] for j in keep_idx])
    m = len(keep_idx)
    if m == 0:
        return float("nan"), float("nan"), float("nan"), 0
    keep = np.asarray(keep_idx)
    boot = np.empty(B); ninf = 0
    for b in range(B):
        idx = keep[rng.integers(0, m, size=m)]
        _, Lb = agg_null([segs_p[j] for j in idx], [span_p[j] for j in idx])
        if not np.isfinite(Lb):
            ninf += 1; Lb = INF_SENTINEL
        boot[b] = Lb
    lo, hi = np.percentile(boot, CI)
    return float(point), float(lo), float(hi), ninf


def main():
    t0 = time.time()
    args = sys.argv[1:]
    bflag = [a for a in args if a.startswith("--B=")]
    B = int(bflag[0].split("=", 1)[1]) if bflag else B_DEFAULT
    args = [a for a in args if not a.startswith("--B=")]
    chroms = args or [f"chr{n}" for n in range(1, 23)]

    print(f"parsing children (chroms={chroms}, gaps={GAPS}) ...", flush=True)
    segs_by_gap, childtot, childspan, nchild = parse_multigap(chroms)

    # F_ROH per child (at SCREEN_GAP) -> reused for both tables
    froh = {}
    for p in POPS:
        if childtot[p] is None:
            continue
        froh[p] = np.array([childtot[p][j] / childspan[p][j] if childspan[p][j] else 0.0
                            for j in range(nchild[p])])

    rng = np.random.default_rng(RNG_SEED)

    # ---- TABLE A: gap x screen surface ----
    surfaceA = []   # (pop, gap, screen_label, n_used, point, lo, hi, ninf)
    for p in POPS:
        if p not in froh:
            continue
        all_idx = list(range(nchild[p]))
        scr_idx = [j for j in all_idx if froh[p][j] <= m21.OUTLIER_F]
        for g in GAPS:
            for label, idx in (("all", all_idx), ("screened", scr_idx)):
                pt, lo, hi, ninf = lstar_ci(segs_by_gap[g][p], childspan[p], idx, rng, B)
                surfaceA.append((p, g, label, len(idx), pt, lo, hi, ninf))
                print(f"  A {p} gap={g} {label:8s} n={len(idx):3d} "
                      f"L*={pt:.2f} CI[{lo:.2f},{hi:.2f}] nonconv={ninf}", flush=True)

    # ---- TABLE B: screen-cutoff sweep at the operating gap ----
    surfaceB = []   # (pop, cutoff_label, n_retained, point, lo, hi, ninf)
    for p in POPS:
        if p not in froh:
            continue
        for cut in OUTLIER_LADDER:
            if cut is None:
                idx = list(range(nchild[p])); clabel = "none"
            else:
                idx = [j for j in range(nchild[p]) if froh[p][j] <= cut]
                clabel = f"{cut:.4f}"
            pt, lo, hi, ninf = lstar_ci(segs_by_gap[OPERATING_GAP][p], childspan[p], idx, rng, B)
            surfaceB.append((p, clabel, len(idx), pt, lo, hi, ninf))
            print(f"  B {p} OUTLIER_F={clabel:7s} n={len(idx):3d} "
                  f"L*={pt:.2f} CI[{lo:.2f},{hi:.2f}]", flush=True)

    with OUT_TSV.open("w", encoding="utf-8") as fh:
        fh.write("# TABLE_A gap x screen surface; bootstrap unit=trio child; "
                 f"B={B} seed={RNG_SEED} CI={CI[0]}-{CI[1]}pct screen_gap={SCREEN_GAP}\n")
        fh.write("table\tpop\tgap_tol\tscreen\tn_children\tLstar_Mb\tCI_lo_Mb\tCI_hi_Mb\tn_nonconverge\n")
        for p, g, lab, n, pt, lo, hi, ninf in surfaceA:
            fh.write(f"A\t{p}\t{g}\t{lab}\t{n}\t{pt:.2f}\t{lo:.2f}\t{hi:.2f}\t{ninf}\n")
        fh.write(f"# TABLE_B screen-cutoff sweep at gap={OPERATING_GAP}; "
                 "n_retained shows how many children survive each F_ROH cutoff\n")
        fh.write("table\tpop\tOUTLIER_F\tscreen\tn_retained\tLstar_Mb\tCI_lo_Mb\tCI_hi_Mb\tn_nonconverge\n")
        for p, clab, n, pt, lo, hi, ninf in surfaceB:
            fh.write(f"B\t{p}\t{clab}\t-\t{n}\t{pt:.2f}\t{lo:.2f}\t{hi:.2f}\t{ninf}\n")

    lines = [
        "# Decisive ROH length L* as a surface (item-2 hardening)",
        f"# chroms={','.join(chroms)} children={sum(nchild.values())} B={B} "
        f"seed={RNG_SEED} CI={CI[0]}-{CI[1]}pct screen_gap={SCREEN_GAP} wall={time.time()-t0:.0f}s",
        f"# OUTLIER_F(screen, Table A)={m21.OUTLIER_F}  prior pi={m21.PI}  posterior>={m21.T_DEC}",
        "",
        "TABLE A -- L* (Mb) [95% CI] over gap tolerance x relatedness screen:",
        "pop      gap   all (n)              screened (n)",
    ]
    bykey = {(p, g, lab): (n, pt, lo, hi) for p, g, lab, n, pt, lo, hi, _ in surfaceA}
    for p in POPS:
        if p not in froh:
            continue
        for g in GAPS:
            na, pa, loa, hia = bykey[(p, g, "all")]
            ns, ps, los, his = bykey[(p, g, "screened")]
            lines.append(f"{p:8s} {g:<4d} {pa:.2f} [{loa:.2f},{hia:.2f}] (n={na:<3d})   "
                         f"{ps:.2f} [{los:.2f},{his:.2f}] (n={ns:<3d})")
    lines += [
        "",
        f"TABLE B -- L* (Mb) [95% CI] and retained n vs screen cutoff, at gap={OPERATING_GAP}:",
        "pop      OUTLIER_F   n_retained   L* [95% CI]",
    ]
    for p, clab, n, pt, lo, hi, ninf in surfaceB:
        lines.append(f"{p:8s} {clab:9s}   {n:<10d}   {pt:.2f} [{lo:.2f}, {hi:.2f}]")
    lines += [
        "",
        "Reading: the gap axis (Table A) shows the headline L* is an operating-point",
        "choice -- gap=0 (no tolerated mismatch) ~ the analytic strict bound, gap=1 the",
        "reported ~1.5 Mb, gap=2 ~2-2.6 Mb -- matching the run-merging tolerance real ROH",
        "callers expose (cf. PLINK --homozyg-gap). The screen axis / Table B show how much",
        "L* and the retained sample size move with the cryptic-relatedness cutoff; where a",
        "large fraction of a population is removed (e.g. AMR), the screened cell rests on a",
        "small self-selected subset and should be read with its n.",
    ]
    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}\n  total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
