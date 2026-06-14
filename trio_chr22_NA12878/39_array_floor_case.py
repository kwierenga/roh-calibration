"""
39_array_floor_case.py - re-analyse the genome-wide platform chunks (script 34
checkpoints in _platform_chunks/) to test, directly from committed data, whether
the calibrated ROH evidence floor transfers from WGS to SNP arrays.

No VCF re-parse: reuses the 22 per-chrom pickles. For each platform x pop it
computes the decisive length L* under three screens:

  L_all     : no relatedness screen
  L_abs     : absolute F_ROH > 0.0156 screen (the WGS-tuned cutoff, as deployed)
  L_qmatch  : drop the SAME NUMBER of highest-F_ROH children as the WGS_full
              absolute screen drops for that population (a density-robust,
              quantile-matched screen) -- isolates "if we could screen as well
              as WGS, what is the array floor?"

Three questions:
  Q1  Ascertainment vs random at matched density: does top-MAF selection
      (CytoSNP-850K-class, 280/Mb) recover the WGS scale relative to the
      Random-300/Mb foil?
  Q2  Density invariance within top-MAF arrays: is L* flat across
      280 -> 200 -> 100 SNP/Mb (above the LD-tagging floor)?
  Q3  Is the array screen failure a true relatedness signal or an F_ROH
      estimation artifact of sparse density? (compare F_ROH inflation and
      L_qmatch vs L_abs)

Usage:  python 39_array_floor_case.py
"""
import importlib.util
import pickle
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
CHUNK_DIR = HERE / "_platform_chunks"

_spec = importlib.util.spec_from_file_location(
    "trionull21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m21)
POPS = m21.POPS
OUTLIER_F = m21.OUTLIER_F

PLATFORMS = ["WGS_full", "thin_1/2", "thin_1/5", "thin_1/10", "thin_1/20",
             "CytoSNP-850K-class", "CytoScan-HD-class", "250K-class",
             "Random-300/Mb"]


def load_accum():
    """Merge all per-chrom chunks -> accum[plat][pop][child] = {segs,span,burden},
    snps_total[plat][pop], span_total[pop]."""
    accum = {pl: {p: {} for p in POPS} for pl in PLATFORMS}
    snps_total = {pl: {p: 0 for p in POPS} for pl in PLATFORMS}
    span_total = {p: 0.0 for p in POPS}
    paths = sorted(CHUNK_DIR.glob("chr*.pkl"),
                   key=lambda x: int(x.stem[3:]))
    for path in paths:
        with open(path, "rb") as fh:
            ck = pickle.load(fh)
        for pl in PLATFORMS:
            for p in POPS:
                snps_total[pl][p] += ck["snps_used"][pl][p]
                for name, b in ck["chunk"][pl][p].items():
                    a = accum[pl][p].setdefault(
                        name, {"segs": [], "span": 0.0, "burden": 0.0})
                    a["segs"].extend(b["segs"])
                    a["span"] += b["span"]
                    a["burden"] += b["burden"]
        for p in POPS:
            span_total[p] += ck["span_seen"][p]
        del ck
    return accum, snps_total, span_total


def lstar(kids_d, drop_idx):
    keys = list(kids_d.keys())
    keep = [i for i in range(len(keys)) if i not in drop_idx]
    segs = [kids_d[keys[i]]["segs"] for i in keep]
    spans = [kids_d[keys[i]]["span"] for i in keep]
    _, L = m21.agg_null(segs, spans)
    return L


def froh_vec(kids_d):
    keys = list(kids_d.keys())
    return keys, np.array([
        (kids_d[k]["burden"] / kids_d[k]["span"]) if kids_d[k]["span"] else 0.0
        for k in keys])


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print("loading 22 chunks ...", flush=True)
    accum, snps_total, span_total = load_accum()

    # WGS_full absolute-screen drop counts per pop -> quantile-matched target
    wgs_drop = {}
    for p in POPS:
        _, f = froh_vec(accum["WGS_full"][p])
        wgs_drop[p] = int((f > OUTLIER_F).sum())

    print(f"\nWGS_full absolute-screen removals (quantile-match target n): "
          f"{wgs_drop}\n")

    hdr = (f"{'platform':<20}{'pop':<5}{'dens':>6}{'n':>5}{'nflag':>6}"
           f"{'medF':>8}{'p90F':>8}{'L_all':>8}{'L_abs':>8}{'L_qm':>8}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for pl in PLATFORMS:
        for p in POPS:
            kd = accum[pl][p]
            if not kd:
                continue
            keys, f = froh_vec(kd)
            n = len(keys)
            dens = snps_total[pl][p] / span_total[p] if span_total[p] else 0
            nflag = int((f > OUTLIER_F).sum())
            medF, p90 = float(np.median(f)), float(np.percentile(f, 90))
            # absolute screen
            drop_abs = set(np.flatnonzero(f > OUTLIER_F).tolist())
            L_abs = lstar(kd, drop_abs)
            # quantile-matched: drop top wgs_drop[p] by F_ROH
            k = min(wgs_drop[p], n)
            drop_qm = set(np.argsort(f)[::-1][:k].tolist()) if k else set()
            L_qm = lstar(kd, drop_qm)
            L_all = lstar(kd, set())
            row = dict(platform=pl, pop=p, dens=dens, n=n, nflag=nflag,
                       medF=medF, p90=p90, L_all=L_all, L_abs=L_abs, L_qm=L_qm)
            rows.append(row)
            def fmt(x):
                return "inf" if x == float("inf") else f"{x:.2f}"
            print(f"{pl:<20}{p:<5}{dens:>6.0f}{n:>5}{nflag:>6}"
                  f"{medF:>8.4f}{p90:>8.4f}{fmt(L_all):>8}"
                  f"{fmt(L_abs):>8}{fmt(L_qm):>8}")
        print()

    # ---- write machine-readable summary
    out = HERE / "array_floor_case.tsv"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("platform\tpop\tdensity\tn\tn_flag_absscreen\tmedian_Froh\t"
                 "p90_Froh\tL_all\tL_abs_screen\tL_qmatch_screen\n")
        for r in rows:
            def f(x):
                return "inf" if x == float("inf") else f"{x:.3f}"
            fh.write(f"{r['platform']}\t{r['pop']}\t{r['dens']:.0f}\t{r['n']}\t"
                     f"{r['nflag']}\t{r['medF']:.4f}\t{r['p90']:.4f}\t"
                     f"{f(r['L_all'])}\t{f(r['L_abs'])}\t{f(r['L_qm'])}\n")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
