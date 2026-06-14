"""
31_platform_genomewide.py - genome-wide platform sensitivity of the decisive
ROH length L* (supersedes the chr22 prototype, script 23).

WHY: Clinical platforms differ in two ways that both affect the calibrated
minimum-callable ROH length -- (1) DENSITY (arrays have ~100-300 SNPs/Mb of
common variants vs ~3000+ for dense WGS; lower density -> fewer mismatches per
Mb -> longer chance-IBS runs -> larger L*); (2) ASCERTAINMENT (array SNPs are
MAF-selected to be highly informative, partially offsetting the density penalty).
Script 23 demonstrated this on chr22 only. This is the genome-wide version,
producing the per-platform L* needed for clinical lookup tables.

PLATFORMS scored in one pass per (pop, child):
  - WGS_full              dense WGS, all common SNPs (MAF >= 0.05)
  - thin_1/K   K=2,5,10,20  random thinning of common SNPs
  - CytoSNP-850K-class    ~280 SNPs/Mb, top-MAF selection
  - CytoScan-HD-class     ~200 SNPs/Mb, top-MAF selection
  - 250K-class            ~100 SNPs/Mb, top-MAF selection (Hildebrandt-era)
  - Random-300/Mb         array-like density without MAF ascertainment (foil)

For each platform we accumulate per-child segments + span + F_ROH burden across
all 22 autosomes, apply the same F_ROH outlier screen as script 21 (per
platform, so children outlying on one platform may not be outliers on another),
and compute L*_screened from the pooled non-outlier children's segments.

Reuses script 21's roh_lengths + agg_null (single source of truth for the
empirical-null aggregator).

Outputs:
  platform_genomewide.tsv  per (platform, pop) L*_screened + n
  platform_genomewide.txt  human-readable summary table

Usage:
  python 31_platform_genomewide.py chr22   # smoke test, fast
  python 31_platform_genomewide.py         # all 22 autosomes (background, ~2-3h)
"""

import gzip
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"

# import script 21's helpers + constants (single source of truth)
_spec = importlib.util.spec_from_file_location(
    "trionull21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m21)

POPS = m21.POPS
MAF_MIN = m21.MAF_MIN
MIN_KEEP_MB = m21.MIN_KEEP_MB
FROH_MIN_MB = m21.FROH_MIN_MB
OUTLIER_F = m21.OUTLIER_F
AF_PRE = m21.AF_PRE

# platform conditions (in addition to full WGS)
THIN_FACTORS = [2, 5, 10, 20]
ARRAY_TARGETS = [   # (label, target SNPs/Mb, selection mode)
    ("CytoSNP-850K-class",  280, "highmaf"),
    ("CytoScan-HD-class",   200, "highmaf"),
    ("250K-class",          100, "highmaf"),  # Hildebrandt 2009 platform
    ("Random-300/Mb",       300, "random"),
]
RNG_SEED = 17
OUT_TSV = HERE / "platform_genomewide.tsv"
OUT_TXT = HERE / "platform_genomewide.txt"


def parse_chrom(chrom, kids):
    """Read one chrom; return per-pop {pos, hom matrix, maf, child_names}."""
    vcf = (HERE / "chr22_phased.vcf.gz" if chrom == "chr22"
           else DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz")
    if not vcf.exists():
        return None
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break
    colpop = {p: [] for p in POPS}
    names = {p: [] for p in POPS}
    for j, s in enumerate(samples):
        if s in kids:
            colpop[kids[s]].append(j)
            names[kids[s]].append(s)
    rows = {p: [] for p in POPS}
    pos = {p: [] for p in POPS}
    maf = {p: [] for p in POPS}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
                continue
            caf = {}
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
                            caf[p] = min(af, 1 - af)
                        break
            if not caf:
                continue
            gts = f[9:]
            p1 = int(f[1])
            for p, m in caf.items():
                rows[p].append(bytes(1 if gts[c][0] == gts[c][2] else 0
                                     for c in colpop[p]))
                pos[p].append(p1)
                maf[p].append(m)
    out = {}
    for p in POPS:
        if not rows[p]:
            continue
        ncol = len(colpop[p])
        mat = np.frombuffer(b"".join(rows[p]),
                            dtype=np.int8).reshape(len(rows[p]), ncol)
        out[p] = {
            "pos": np.asarray(pos[p], dtype=np.int64),
            "hom": mat,
            "maf": np.asarray(maf[p]),
            "names": names[p],
        }
    return out


def select_indices(n_snps, span_mb, mafs, density, mode, rng):
    """Pick SNP indices for a target density, by random or top-MAF selection."""
    n_target = min(int(density * span_mb), n_snps)
    if n_target <= 0:
        return np.empty(0, dtype=np.int64)
    if mode == "random":
        idx = rng.choice(n_snps, n_target, replace=False)
    else:  # "highmaf"
        idx = np.argsort(mafs)[::-1][:n_target]
    return np.sort(idx)


def main():
    t0 = time.time()
    chroms = sys.argv[1:] or [f"chr{n}" for n in range(1, 23)]

    sp = m21.load_superpop()
    kids = m21.load_children(sp)
    rng = np.random.default_rng(RNG_SEED)

    platforms = ["WGS_full"]
    for k in THIN_FACTORS:
        platforms.append(f"thin_1/{k}")
    for label, _, _ in ARRAY_TARGETS:
        platforms.append(label)

    # accumulator: accum[platform][pop][child_name] = {segs, span, burden}
    accum = {plat: {p: {} for p in POPS} for plat in platforms}
    # for reporting achieved density
    snps_used = {plat: {p: 0.0 for p in POPS} for plat in platforms}
    span_seen = {p: 0.0 for p in POPS}

    for chrom in chroms:
        ts = time.time()
        data = parse_chrom(chrom, kids)
        if data is None:
            print(f"  [{chrom}] SKIP", flush=True)
            continue
        for p in POPS:
            if p not in data:
                continue
            d = data[p]
            pa, hom, mafa, names = d["pos"], d["hom"], d["maf"], d["names"]
            n_snps = pa.size
            ncol = hom.shape[1]
            span = (pa[-1] - pa[0]) / 1e6
            span_seen[p] += span

            # per-platform index sets for this chrom
            sel = {"WGS_full": np.arange(n_snps)}
            for k in THIN_FACTORS:
                sel[f"thin_1/{k}"] = np.arange(0, n_snps, k)
            for label, dens, mode in ARRAY_TARGETS:
                sel[label] = select_indices(n_snps, span, mafa, dens, mode, rng)

            for plat, idx in sel.items():
                if idx.size == 0:
                    continue
                snps_used[plat][p] += idx.size
                pa_sel = pa[idx]
                hom_sel = hom[idx, :]
                for j in range(ncol):
                    name = names[j]
                    bucket = accum[plat][p].setdefault(
                        name, {"segs": [], "span": 0.0, "burden": 0.0})
                    segs = m21.roh_lengths(hom_sel[:, j].astype(bool), pa_sel)
                    bucket["burden"] += float(segs[segs >= FROH_MIN_MB].sum())
                    keep = segs[segs > MIN_KEEP_MB]
                    if keep.size:
                        bucket["segs"].append(keep)
                    bucket["span"] += span
        msg = " ".join(
            f"{p}:{len(data[p]['names']) if p in data else 0}" for p in POPS)
        print(f"  [{chrom}] {msg}  ({time.time()-ts:.0f}s, "
              f"{time.time()-t0:.0f}s total)", flush=True)

    # compute L* per (platform, pop)
    results = {}
    for plat in platforms:
        for p in POPS:
            kids_d = accum[plat][p]
            if not kids_d:
                continue
            n = len(kids_d)
            keys = list(kids_d.keys())
            f_roh = np.array([
                (kids_d[k]["burden"] / kids_d[k]["span"])
                if kids_d[k]["span"] else 0.0 for k in keys])
            out_mask = f_roh > OUTLIER_F
            segs_all = [kids_d[k]["segs"] for k in keys]
            spans_all = [kids_d[k]["span"] for k in keys]
            _, L_all = m21.agg_null(segs_all, spans_all)
            keep = [i for i in range(n) if not out_mask[i]]
            _, L_scr = m21.agg_null([segs_all[i] for i in keep],
                                    [spans_all[i] for i in keep])
            density = (snps_used[plat][p] / span_seen[p]
                       if span_seen[p] else 0.0)
            results[(plat, p)] = {
                "n_all": n, "n_out": int(out_mask.sum()),
                "n_scr": int(n - out_mask.sum()),
                "L_all": L_all, "L_scr": L_scr, "density": density,
            }

    # outputs
    lines = [
        "# Genome-wide platform sensitivity of L* (supersedes chr22 prototype, script 23)",
        f"# chroms={','.join(chroms)}  PI={m21.PI}  GAP_TOL={m21.GAP_TOL}  "
        f"OUTLIER_F={OUTLIER_F}",
        "# L*_screened: empirical decisive ROH length after F_ROH outlier screen "
        "(per platform).",
        f"# wall={time.time()-t0:.0f}s\n",
        "platform\t" + "\t".join(POPS) + "\tavg density (SNPs/Mb)\tnote",
    ]
    for plat in platforms:
        cells, dens_p = [], []
        for p in POPS:
            r = results.get((plat, p))
            cells.append(f"{r['L_scr']:.2f}" if r else "NA")
            if r:
                dens_p.append(r["density"])
        avg_dens = np.mean(dens_p) if dens_p else 0.0
        if plat == "WGS_full":
            note = "dense WGS, full common-SNP density"
        elif plat.startswith("thin"):
            k = int(plat.split("/")[1])
            note = f"random thinning, every {k}th common SNP"
        elif plat == "Random-300/Mb":
            note = "array density without MAF ascertainment (foil)"
        else:
            note = "clinical array proxy, top-MAF ascertainment"
        lines.append(f"{plat}\t" + "\t".join(cells) +
                     f"\t{avg_dens:.0f}\t{note}")

    tsv = ["platform\tpop\tn_children\tn_outliers\tn_screened\t"
           "L_screened_Mb\tL_all_Mb\tdensity_SNPs_per_Mb"]
    for plat in platforms:
        for p in POPS:
            r = results.get((plat, p))
            if not r:
                continue
            tsv.append(f"{plat}\t{p}\t{r['n_all']}\t{r['n_out']}\t{r['n_scr']}\t"
                       f"{r['L_scr']:.2f}\t{r['L_all']:.2f}\t{r['density']:.0f}")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    OUT_TSV.write_text("\n".join(tsv) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}")


if __name__ == "__main__":
    main()
