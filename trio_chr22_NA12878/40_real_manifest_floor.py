"""
40_real_manifest_floor.py - test the SNP-array evidence floor using REAL array
manifests (Wrayner b38 strand files) intersected with 1000G WGS, instead of the
top-MAF proxy of script 34.

For each real array we select the 1000G common SNPs (MAF>=0.05) whose GRCh38
position is an actual probe on that array, call ROH in the trio children, and
compute the decisive length L* under three screens (no / absolute F_ROH>0.0156 /
rank-matched), exactly as script 39. We also score WGS_full and a top-MAF proxy
matched to each array's realized common-SNP count, so real-vs-proxy is apples to
apples.

Build check: 1000G is GRCh38; if a strand file were secretly hg19 the overlap
with common-SNP positions would collapse. We report overlap % per array.

Usage:
  python 40_real_manifest_floor.py chr22            # smoke test (one chrom)
  python 40_real_manifest_floor.py                  # all 22, per-chrom checkpoint
  python 40_real_manifest_floor.py --aggregate-only
"""
import gc
import importlib.util
import pickle
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
STRAND_DIR = HERE / "external" / "array_strand"
CHUNK_DIR = HERE / "_manifest_chunks"

_s21 = importlib.util.spec_from_file_location("trionull21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_s21); _s21.loader.exec_module(m21)
_s34 = importlib.util.spec_from_file_location("platform34", HERE / "34_platform_genomewide_v2.py")
m34 = importlib.util.module_from_spec(_s34); _s34.loader.exec_module(m34)

POPS = m21.POPS
OUTLIER_F = m21.OUTLIER_F

# real arrays: label -> strand filename
ARRAYS = {
    "CytoSNP-850K": "CytoSNP-850K_v1-1_iScan_A1-b38.strand",
    "GSA-24v3":     "GSA-24v3-0_A2-b38.strand",
    "MEGA":         "Multi-EthnicGlobal_A1-b38.strand",
}
RNG_SEED = 17


def load_manifest(fname):
    """strand file cols: SNPid chr pos %match strand TOP. Return {chrom_int: sorted np.array(pos)}."""
    by_chrom = {c: [] for c in range(1, 23)}
    path = STRAND_DIR / fname
    with open(path, "r") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3:
                continue
            c = f[1]
            if not c.isdigit():
                continue
            ci = int(c)
            if 1 <= ci <= 22:
                try:
                    by_chrom[ci].append(int(f[2]))
                except ValueError:
                    continue
    return {c: np.array(sorted(set(v)), dtype=np.int64) for c, v in by_chrom.items() if v}


def process_chrom(chrom, kids, manifests, rng):
    data = m34.parse_chrom(chrom, kids)
    if data is None:
        return None
    ci = int(chrom.replace("chr", ""))
    platforms = ["WGS_full"] + list(ARRAYS) + [a + "_proxy" for a in ARRAYS]
    chunk = {pl: {p: {} for p in POPS} for pl in platforms}
    snps_used = {pl: {p: 0 for p in POPS} for pl in platforms}
    overlap = {a: {"array_on_chrom": 0, "common_hit": 0} for a in ARRAYS}
    span_seen = {p: 0.0 for p in POPS}

    for p in POPS:
        if p not in data:
            continue
        d = data[p]
        pa, hom, mafa, names = d["pos"], d["hom"], d["maf"], d["names"]
        n_snps = pa.size
        span = (pa[-1] - pa[0]) / 1e6
        span_seen[p] = span

        sel = {"WGS_full": np.arange(n_snps)}
        for a, _ in ARRAYS.items():
            mpos = manifests[a].get(ci)
            if mpos is None or mpos.size == 0:
                sel[a] = np.empty(0, dtype=np.int64)
                sel[a + "_proxy"] = np.empty(0, dtype=np.int64)
                continue
            mask = np.isin(pa, mpos)               # common SNP positions that are on the array
            idx = np.flatnonzero(mask)
            sel[a] = idx
            if p == POPS[0]:
                overlap[a]["array_on_chrom"] += int(mpos.size)
                overlap[a]["common_hit"] += int(idx.size)
            # proxy: top-MAF at the SAME realized count
            ntop = idx.size
            sel[a + "_proxy"] = (np.sort(np.argsort(mafa)[::-1][:ntop])
                                 if ntop else np.empty(0, dtype=np.int64))

        for pl, idx in sel.items():
            if idx.size == 0:
                continue
            snps_used[pl][p] = int(idx.size)
            pa_sel, hom_sel = pa[idx], hom[idx, :]
            for j in range(hom.shape[1]):
                segs = m21.roh_lengths(hom_sel[:, j].astype(bool), pa_sel)
                burden = float(segs[segs >= m21.FROH_MIN_MB].sum())
                keep = segs[segs > m21.MIN_KEEP_MB]
                b = chunk[pl][p].setdefault(names[j], {"segs": [], "span": 0.0, "burden": 0.0})
                if keep.size:
                    b["segs"].append(keep.astype(np.float32))
                b["span"] += span
                b["burden"] += burden
            del hom_sel
        del data[p]["hom"]
    del data
    gc.collect()
    return {"chunk": chunk, "snps_used": snps_used, "overlap": overlap,
            "span_seen": span_seen, "platforms": platforms}


def lstar(kids_d, drop_idx):
    keys = list(kids_d.keys())
    keep = [i for i in range(len(keys)) if i not in drop_idx]
    segs = [kids_d[keys[i]]["segs"] for i in keep]
    spans = [kids_d[keys[i]]["span"] for i in keep]
    return m21.agg_null(segs, spans)[1]


def aggregate_and_report(chunk_paths):
    first = pickle.load(open(chunk_paths[0], "rb"))
    platforms = first["platforms"]
    accum = {pl: {p: {} for p in POPS} for pl in platforms}
    snps_total = {pl: {p: 0 for p in POPS} for pl in platforms}
    span_total = {p: 0.0 for p in POPS}
    ov = {a: {"array_on_chrom": 0, "common_hit": 0} for a in ARRAYS}
    for path in chunk_paths:
        ck = pickle.load(open(path, "rb"))
        for pl in platforms:
            for p in POPS:
                snps_total[pl][p] += ck["snps_used"][pl][p]
                for name, b in ck["chunk"][pl][p].items():
                    a = accum[pl][p].setdefault(name, {"segs": [], "span": 0.0, "burden": 0.0})
                    a["segs"].extend(b["segs"]); a["span"] += b["span"]; a["burden"] += b["burden"]
        for p in POPS:
            span_total[p] += ck["span_seen"][p]
        for a in ARRAYS:
            ov[a]["array_on_chrom"] += ck["overlap"][a]["array_on_chrom"]
            ov[a]["common_hit"] += ck["overlap"][a]["common_hit"]
        del ck; gc.collect()

    # WGS drop counts -> quantile-match target
    wgs_drop = {}
    for p in POPS:
        kd = accum["WGS_full"][p]
        f = np.array([(kd[k]["burden"]/kd[k]["span"]) if kd[k]["span"] else 0.0 for k in kd])
        wgs_drop[p] = int((f > OUTLIER_F).sum())

    print("\n# Real-array overlap with 1000G common SNPs (build check; low overlap => wrong build):")
    for a in ARRAYS:
        print(f"  {a}: {ov[a]['common_hit']} common-SNP hits of {ov[a]['array_on_chrom']} "
              f"autosomal probes ({100*ov[a]['common_hit']/max(1,ov[a]['array_on_chrom']):.1f}%)")

    hdr = f"\n{'platform':<20}{'pop':<5}{'dens':>6}{'n':>5}{'L_all':>8}{'L_abs':>8}{'L_qm':>8}"
    print(hdr); print("-" * len(hdr))
    rows = []
    for pl in platforms:
        for p in POPS:
            kd = accum[pl][p]
            if not kd:
                continue
            keys = list(kd.keys()); n = len(keys)
            f = np.array([(kd[k]["burden"]/kd[k]["span"]) if kd[k]["span"] else 0.0 for k in keys])
            dens = snps_total[pl][p]/span_total[p] if span_total[p] else 0
            drop_abs = set(np.flatnonzero(f > OUTLIER_F).tolist())
            k = min(wgs_drop[p], n)
            drop_qm = set(np.argsort(f)[::-1][:k].tolist()) if k else set()
            L_all, L_abs, L_qm = lstar(kd, set()), lstar(kd, drop_abs), lstar(kd, drop_qm)
            rows.append((pl, p, dens, n, L_all, L_abs, L_qm))
            fmt = lambda x: "inf" if x == float("inf") else f"{x:.2f}"
            print(f"{pl:<20}{p:<5}{dens:>6.0f}{n:>5}{fmt(L_all):>8}{fmt(L_abs):>8}{fmt(L_qm):>8}")
        print()

    out = HERE / "real_manifest_floor.tsv"
    with open(out, "w") as fh:
        fh.write("platform\tpop\tdensity\tn\tL_all\tL_abs_screen\tL_qmatch_screen\n")
        for pl, p, dens, n, a, ab, qm in rows:
            g = lambda x: "inf" if x == float("inf") else f"{x:.3f}"
            fh.write(f"{pl}\t{p}\t{dens:.0f}\t{n}\t{g(a)}\t{g(ab)}\t{g(qm)}\n")
    print(f"  -> {out}")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    t0 = time.time()
    args = [a for a in sys.argv[1:] if a != "--aggregate-only"]
    agg_only = "--aggregate-only" in sys.argv
    chroms = args or [f"chr{n}" for n in range(1, 23)]
    CHUNK_DIR.mkdir(exist_ok=True)

    if not agg_only:
        print("loading manifests ...", flush=True)
        manifests = {a: load_manifest(f) for a, f in ARRAYS.items()}
        for a in ARRAYS:
            tot = sum(v.size for v in manifests[a].values())
            print(f"  {a}: {tot} autosomal probe positions", flush=True)
        sp = m21.load_superpop(); kids = m21.load_children(sp)
        rng = np.random.default_rng(RNG_SEED)
        for chrom in chroms:
            cp = CHUNK_DIR / f"{chrom}.pkl"
            if cp.exists():
                print(f"  [{chrom}] chunk exists, skip", flush=True); continue
            ts = time.time()
            r = process_chrom(chrom, kids, manifests, rng)
            if r is None:
                print(f"  [{chrom}] SKIP (no VCF)", flush=True); continue
            pickle.dump(r, open(cp, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
            print(f"  [{chrom}] done ({time.time()-ts:.0f}s, {time.time()-t0:.0f}s total)", flush=True)
            del r; gc.collect()

    chunk_paths = sorted(CHUNK_DIR.glob("chr*.pkl"),
                         key=lambda x: int(x.stem.replace("chr", "")))
    if chunk_paths:
        print(f"\naggregating {len(chunk_paths)} chunk(s) ...", flush=True)
        aggregate_and_report(chunk_paths)


if __name__ == "__main__":
    main()
