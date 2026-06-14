"""
34_platform_genomewide_v2.py - genome-wide platform sensitivity of L*, with
per-chromosome on-disk checkpointing.

Supersedes script 31. Script 31 finished chr1 cleanly (1030s) and was then
OS-killed at the start of chr2 with no Python exception in the log -- the
in-memory accumulator (per-platform per-pop per-child list of small numpy
segment arrays) plus chr1's transient hom_sel copies hit the working-set limit.

Fix: at the end of each chromosome's processing, pickle the in-memory
accumulator to `_platform_chunks/chr{N}.pkl` and reset the accumulator to
empty. After all chroms, load the 22 chunks one at a time and merge them
into a single aggregate, then compute L* per (platform, pop). Peak memory is
bounded to "one chrom's data" rather than "cumulative all-chroms data".

Same platform set as script 31:
  - WGS_full
  - thin_1/{2, 5, 10, 20}
  - CytoSNP-850K-class (~280 SNPs/Mb, top-MAF)
  - CytoScan-HD-class  (~200 SNPs/Mb, top-MAF)
  - 250K-class         (~100 SNPs/Mb, top-MAF) -- Hildebrandt-era
  - Random-300/Mb      (no MAF ascertainment, foil)

Outputs:
  platform_genomewide.tsv  per (platform, pop) L*_screened + density
  platform_genomewide.txt  human-readable summary
  _platform_chunks/        per-chrom pickle checkpoints (kept; small)

Usage:
  python 34_platform_genomewide_v2.py chr22     # smoke test
  python 34_platform_genomewide_v2.py           # all 22 autosomes, background
  python 34_platform_genomewide_v2.py --aggregate-only   # skip parse, just aggregate existing chunks
"""

import gc
import gzip
import importlib.util
import pickle
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"
CHUNK_DIR = HERE / "_platform_chunks"

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

THIN_FACTORS = [2, 5, 10, 20]
ARRAY_TARGETS = [
    ("CytoSNP-850K-class",  280, "highmaf"),
    ("CytoScan-HD-class",   200, "highmaf"),
    ("250K-class",          100, "highmaf"),
    ("Random-300/Mb",       300, "random"),
]
RNG_SEED = 17
OUT_TSV = HERE / "platform_genomewide.tsv"
OUT_TXT = HERE / "platform_genomewide.txt"


def parse_chrom(chrom, kids):
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
    n_target = min(int(density * span_mb), n_snps)
    if n_target <= 0:
        return np.empty(0, dtype=np.int64)
    if mode == "random":
        idx = rng.choice(n_snps, n_target, replace=False)
    else:
        idx = np.argsort(mafs)[::-1][:n_target]
    return np.sort(idx)


def process_chrom(chrom, kids, rng, platforms):
    """Return per-chrom accum_chunk: {plat: {pop: {child: dict}}} for one chrom only."""
    data = parse_chrom(chrom, kids)
    if data is None:
        return None
    chunk = {plat: {p: {} for p in POPS} for plat in platforms}
    snps_used = {plat: {p: 0 for p in POPS} for plat in platforms}
    span_seen = {p: 0.0 for p in POPS}

    for p in POPS:
        if p not in data:
            continue
        d = data[p]
        pa, hom, mafa, names = d["pos"], d["hom"], d["maf"], d["names"]
        n_snps = pa.size
        ncol = hom.shape[1]
        span = (pa[-1] - pa[0]) / 1e6
        span_seen[p] = span

        sel = {"WGS_full": np.arange(n_snps)}
        for k in THIN_FACTORS:
            sel[f"thin_1/{k}"] = np.arange(0, n_snps, k)
        for label, dens, mode in ARRAY_TARGETS:
            sel[label] = select_indices(n_snps, span, mafa, dens, mode, rng)

        for plat, idx in sel.items():
            if idx.size == 0:
                continue
            snps_used[plat][p] = int(idx.size)
            pa_sel = pa[idx]
            hom_sel = hom[idx, :]
            for j in range(ncol):
                name = names[j]
                segs = m21.roh_lengths(hom_sel[:, j].astype(bool), pa_sel)
                burden = float(segs[segs >= FROH_MIN_MB].sum())
                keep = segs[segs > MIN_KEEP_MB]
                bucket = chunk[plat][p].setdefault(
                    name, {"segs": [], "span": 0.0, "burden": 0.0})
                if keep.size:
                    bucket["segs"].append(keep.astype(np.float32))
                bucket["span"] += span
                bucket["burden"] += burden
            del hom_sel
        del data[p]["hom"]
    del data
    gc.collect()
    return {"chunk": chunk, "snps_used": snps_used, "span_seen": span_seen}


def aggregate(platforms, chunk_paths):
    """Merge per-chrom pickles into the final accumulator."""
    accum = {plat: {p: {} for p in POPS} for plat in platforms}
    snps_total = {plat: {p: 0 for p in POPS} for plat in platforms}
    span_total = {p: 0.0 for p in POPS}
    for path in chunk_paths:
        with open(path, "rb") as fh:
            ck = pickle.load(fh)
        for plat in platforms:
            for p in POPS:
                snps_total[plat][p] += ck["snps_used"][plat][p]
                for name, bucket in ck["chunk"][plat][p].items():
                    a = accum[plat][p].setdefault(
                        name, {"segs": [], "span": 0.0, "burden": 0.0})
                    a["segs"].extend(bucket["segs"])
                    a["span"] += bucket["span"]
                    a["burden"] += bucket["burden"]
        for p in POPS:
            span_total[p] += ck["span_seen"][p]
        del ck
        gc.collect()
    return accum, snps_total, span_total


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    t0 = time.time()
    args = sys.argv[1:]
    aggregate_only = "--aggregate-only" in args
    args = [a for a in args if a != "--aggregate-only"]
    chroms = args or [f"chr{n}" for n in range(1, 23)]

    platforms = ["WGS_full"]
    for k in THIN_FACTORS:
        platforms.append(f"thin_1/{k}")
    for label, _, _ in ARRAY_TARGETS:
        platforms.append(label)

    CHUNK_DIR.mkdir(exist_ok=True)
    chunk_paths = []
    if not aggregate_only:
        sp = m21.load_superpop()
        kids = m21.load_children(sp)
        rng = np.random.default_rng(RNG_SEED)
        for chrom in chroms:
            chunk_path = CHUNK_DIR / f"{chrom}.pkl"
            if chunk_path.exists():
                print(f"  [{chrom}] CHUNK EXISTS, skip", flush=True)
                chunk_paths.append(chunk_path)
                continue
            ts = time.time()
            r = process_chrom(chrom, kids, rng, platforms)
            if r is None:
                print(f"  [{chrom}] SKIP (no VCF)", flush=True)
                continue
            with open(chunk_path, "wb") as fh:
                pickle.dump(r, fh, protocol=pickle.HIGHEST_PROTOCOL)
            chunk_paths.append(chunk_path)
            msg = " ".join(
                f"{p}:{len(r['chunk'][platforms[0]][p])}" for p in POPS)
            print(f"  [{chrom}] {msg}  ({time.time()-ts:.0f}s, "
                  f"{time.time()-t0:.0f}s total, chunk={chunk_path.stat().st_size//1024}KB)",
                  flush=True)
            del r
            gc.collect()
    else:
        chunk_paths = sorted(CHUNK_DIR.glob("chr*.pkl"))
        print(f"aggregate-only mode: {len(chunk_paths)} chunks", flush=True)

    if not chunk_paths:
        print("no chunks to aggregate")
        return

    print("\naggregating chunks ...", flush=True)
    accum, snps_total, span_total = aggregate(platforms, chunk_paths)

    # compute L* per (platform, pop), with F_ROH outlier screen per platform
    results = {}
    for plat in platforms:
        for p in POPS:
            kids_d = accum[plat][p]
            if not kids_d:
                continue
            keys = list(kids_d.keys())
            n = len(keys)
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
            density = (snps_total[plat][p] / span_total[p]
                       if span_total[p] else 0.0)
            results[(plat, p)] = {
                "n_all": n, "n_out": int(out_mask.sum()),
                "n_scr": int(n - out_mask.sum()),
                "L_all": L_all, "L_scr": L_scr, "density": density,
            }

    lines = [
        "# Genome-wide platform sensitivity of L* (v2, per-chrom checkpointed)",
        f"# chroms={','.join(chroms)}  PI={m21.PI}  GAP_TOL={m21.GAP_TOL}  "
        f"OUTLIER_F={OUTLIER_F}",
        f"# wall={time.time()-t0:.0f}s\n",
        "platform\t" + "\t".join(POPS) + "\tavg density\tnote",
    ]
    for plat in platforms:
        cells, dens_p = [], []
        for p in POPS:
            r = results.get((plat, p))
            cells.append(f"{r['L_scr']:.2f}" if r else "NA")
            if r:
                dens_p.append(r["density"])
        avg_dens = float(np.mean(dens_p)) if dens_p else 0.0
        if plat == "WGS_full":
            note = "dense WGS, full common-SNP density"
        elif plat.startswith("thin"):
            k = int(plat.split("/")[1])
            note = f"random thinning, every {k}th common SNP"
        elif plat == "Random-300/Mb":
            note = "array density without MAF ascertainment (foil)"
        else:
            note = "clinical array proxy, top-MAF ascertainment"
        lines.append(f"{plat}\t" + "\t".join(cells)
                     + f"\t{avg_dens:.0f}\t{note}")

    tsv = ["platform\tpop\tn_children\tn_outliers\tn_screened\t"
           "L_screened_Mb\tL_all_Mb\tdensity_SNPs_per_Mb"]
    for plat in platforms:
        for p in POPS:
            r = results.get((plat, p))
            if not r:
                continue
            tsv.append(f"{plat}\t{p}\t{r['n_all']}\t{r['n_out']}\t"
                       f"{r['n_scr']}\t{r['L_scr']:.2f}\t{r['L_all']:.2f}\t"
                       f"{r['density']:.0f}")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    OUT_TSV.write_text("\n".join(tsv) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n  -> {OUT_TSV}\n  -> {OUT_TXT}")


if __name__ == "__main__":
    main()
