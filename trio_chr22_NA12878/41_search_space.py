"""
41 — Search-space cost of lowering the ROH floor (reviewer point 1).

Per-tract FDR control does NOT bound per-CASE false-positive burden. The
operational reason labs use a 5-10 Mb floor is that lowering it multiplies the
number of candidate ROH (and candidate recessive genes) the analyst must triage
per case. This script quantifies that growth directly on the genome-wide
WGS trio-child panels already on disk (_platform_chunks/*.pkl, WGS_full),
the same children the calibration uses as its outbred background.

For each superpopulation and pooled, at floors {1.6, 5, 10} Mb, we report the
per-case number of ROH and total Mb under ROH, for:
  (a) the screened-baseline set (F_ROH <= OUTLIER_F, the calibration background)
  (b) the full set (includes elevated-background individuals) — shows how the
      search space grows with the patient's own autozygosity (reviewer point 7).

Candidate recessive genes are estimated as (Mb under ROH) x PROT_CODING_PER_MB,
a genome-average protein-coding density stated explicitly; the exact, verifiable
quantities are the ROH count and the Mb under ROH.

Pure stdlib + numpy; reuses module 21 constants (OUTLIER_F, FROH_MIN_MB).
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
FROH_MIN_MB = m21.FROH_MIN_MB      # 1.0 Mb: ROH >= this count toward F_ROH burden
OUTLIER_F = m21.OUTLIER_F          # 0.0156 (~2nd-cousin): screen for recent shared ancestry

PLATFORM = "WGS_full"
FLOORS = [1.6, 5.0, 10.0]
# ~20,000 protein-coding genes over ~2,875 Mb of autosome ~= 7.0/Mb (GENCODE-scale,
# stated as a transparent average; gene density is non-uniform so this is an estimate).
PROT_CODING_PER_MB = 7.0

OUT_TSV = HERE / "search_space.tsv"
OUT_TXT = HERE / "search_space.txt"


def aggregate_wgs():
    """Merge per-chrom pickles -> per-(pop, child) genome-wide ROH segment list (Mb)."""
    accum = {p: {} for p in POPS}          # pop -> child -> {"segs":[float], "span":float}
    chunks = sorted(CHUNK_DIR.glob("*.pkl"))
    if not chunks:
        raise SystemExit(f"no panel chunks in {CHUNK_DIR}")
    for cp in chunks:
        with cp.open("rb") as fh:
            ck = pickle.load(fh)
        ck = ck.get("chunk", ck)        # chunk files wrap payload under "chunk"
        plat = ck.get(PLATFORM)
        if plat is None:
            continue
        for p in POPS:
            for child, bucket in plat.get(p, {}).items():
                a = accum[p].setdefault(child, {"segs": [], "span": 0.0})
                for arr in bucket["segs"]:
                    a["segs"].extend(float(x) for x in np.asarray(arr).ravel())
                a["span"] += float(bucket["span"])
    return accum


def summarize(seg_lists):
    """seg_lists: list (one per child) of np arrays of ROH lengths in Mb.
    Returns dict floor -> per-case count array and per-case Mb array."""
    out = {}
    for T in FLOORS:
        counts = np.array([int((s >= T).sum()) for s in seg_lists], dtype=float)
        mb = np.array([float(s[s >= T].sum()) for s in seg_lists], dtype=float)
        out[T] = (counts, mb)
    return out


def q(a, p):
    return float(np.percentile(a, p)) if len(a) else float("nan")


def amax(a):
    return float(np.max(a)) if len(a) else float("nan")


def amean(a):
    return float(np.mean(a)) if len(a) else float("nan")


def fmt_block(label, seg_lists):
    n = len(seg_lists)
    rows = []
    for T in FLOORS:
        counts, mb = summarize(seg_lists)[T]
        genes = mb * PROT_CODING_PER_MB
        rows.append((label, n, T,
                     q(counts, 50), q(counts, 25), q(counts, 75), amean(counts), amax(counts),
                     q(mb, 50), q(mb, 75),
                     q(genes, 50), q(genes, 75)))
    return rows


def main():
    accum = aggregate_wgs()
    for p in POPS:
        print(f"  [load] {p}: {len(accum[p])} children")

    header = ["set", "pop", "n_cases", "floor_Mb",
              "roh_median", "roh_q25", "roh_q75", "roh_mean", "roh_max",
              "mb_median", "mb_q75", "genes_median", "genes_q75"]
    tsv_rows = []
    txt = []
    txt.append("Search-space cost of lowering the ROH floor (WGS trio children)")
    txt.append(f"Platform={PLATFORM}  screen: F_ROH<= {OUTLIER_F} (~2nd-cousin)  "
               f"genes ~= Mb x {PROT_CODING_PER_MB}/Mb (genome-average protein-coding)")
    txt.append("=" * 78)

    # Build per-pop screened / all child segment lists, plus pooled.
    pooled_scr, pooled_all = [], []
    pop_lists = {}
    for p in POPS:
        scr, allc = [], []
        for child, a in accum[p].items():
            segs = np.array(a["segs"], dtype=float)
            span = a["span"] if a["span"] > 0 else 1.0
            froh = float(segs[segs >= FROH_MIN_MB].sum()) / span
            allc.append(segs)
            if froh <= OUTLIER_F:
                scr.append(segs)
        pop_lists[p] = (scr, allc)
        pooled_scr.extend(scr)
        pooled_all.extend(allc)

    def emit(set_label, pop_label, seg_lists):
        for r in fmt_block(pop_label, seg_lists):
            (_, n, T, med, q25, q75, mean, mx, mbmed, mbq75, gmed, gq75) = r
            tsv_rows.append([set_label, pop_label, n, f"{T:.1f}",
                             f"{med:.0f}", f"{q25:.0f}", f"{q75:.0f}", f"{mean:.2f}", f"{mx:.0f}",
                             f"{mbmed:.2f}", f"{mbq75:.2f}", f"{gmed:.0f}", f"{gq75:.0f}"])

    for set_label, getter in (("screened", 0), ("all", 1)):
        txt.append("")
        txt.append(f"[{set_label.upper()} set]")
        txt.append(f"{'pop':>7} {'n':>4} {'floor':>6} {'ROH/case median(IQR)':>22} "
                   f"{'Mb/case med':>11} {'cand.genes med':>14}")
        for p in POPS:
            seg_lists = pop_lists[p][getter]
            emit(set_label, p, seg_lists)
            for T in FLOORS:
                counts, mb = summarize(seg_lists)[T]
                genes = mb * PROT_CODING_PER_MB
                txt.append(f"{p:>7} {len(seg_lists):>4} {T:>6.1f} "
                           f"{q(counts,50):>6.0f} ({q(counts,25):>3.0f}-{q(counts,75):<3.0f})       "
                           f"{q(mb,50):>11.2f} {q(genes,50):>14.0f}")
        pooled = pooled_scr if set_label == "screened" else pooled_all
        emit(set_label, "POOLED", pooled)
        for T in FLOORS:
            counts, mb = summarize(pooled)[T]
            genes = mb * PROT_CODING_PER_MB
            txt.append(f"{'POOLED':>7} {len(pooled):>4} {T:>6.1f} "
                       f"{q(counts,50):>6.0f} ({q(counts,25):>3.0f}-{q(counts,75):<3.0f})       "
                       f"{q(mb,50):>11.2f} {q(genes,50):>14.0f}")

    # Fold growth factor (1.6 vs 10) for the pooled screened set.
    cs = summarize(pooled_scr)
    c16 = cs[1.6][0]; c10 = cs[10.0][0]
    g16 = float(np.median(c16)); g10 = float(np.median(c10))
    g5 = float(np.median(cs[5.0][0]))
    ratio = f"{g16/g10:.1f}x" if g10 else "from ~0"
    txt.append("")
    txt.append(f"Pooled screened (outbred baseline): median ROH/case "
               f"10 Mb={g10:.0f}, 5 Mb={g5:.0f}, 1.6 Mb={g16:.0f} -> "
               f"lowering the floor surfaces ~{g16:.0f} tracts/case that 5-10 Mb shows none of ({ratio}).")

    OUT_TSV.write_text("\t".join(header) + "\n" +
                       "\n".join("\t".join(map(str, r)) for r in tsv_rows) + "\n")
    OUT_TXT.write_text("\n".join(txt) + "\n")
    print("\n".join(txt))
    print(f"\nwrote {OUT_TSV.name}, {OUT_TXT.name}")


if __name__ == "__main__":
    main()
