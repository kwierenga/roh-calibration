"""
44 — Per-locus EMPIRICAL validation (reviewer point 8).

The central thesis — a fixed length is not a fixed weight of evidence, because the
decisive length L* falls as local recombination rises — has so far been shown only
from the closed-form law (H-bar^n_eff), which we elsewhere show is anti-conservative
by ~2-3x. This script tests it EMPIRICALLY: it bins the autosome by local deCODE
recombination rate, measures the real background-run survival p_background(L) within
each bin directly from the cryptic-relatedness-screened trio children, and reports
the empirical decisive length L*(rate) — with a children-bootstrap for the thin-tail
noise. If L* falls monotonically with rate on real data, the locus thesis is
demonstrated, not asserted.

p_background_bin(L) = Sum_seg max(0, len-L) / (n_children * footprint_Mb[bin]);
L*_bin = smallest L with p_background_bin(L) <= THR_PC (posterior>=0.95, 1st-cousin).
Segments are assigned to a rate bin by their own mean cM/Mb (= the n_eff-relevant
rate). Ancestry is pooled: the recombination landscape is population-independent, so
pooling maximizes power for the locus effect.

Usage:  python 44_per_locus_empirical.py chr22     # smoke test
        python 44_per_locus_empirical.py            # all autosomes (background)
"""
import gzip
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"

_s21 = importlib.util.spec_from_file_location("m21", HERE / "21_trio_background_null.py")
m21 = importlib.util.module_from_spec(_s21); _s21.loader.exec_module(m21)
_s16 = importlib.util.spec_from_file_location("m16", HERE / "16_haplotype_ibs_noise.py")
m16 = importlib.util.module_from_spec(_s16); _s16.loader.exec_module(m16)

POPS = m21.POPS
MAF_MIN = m21.MAF_MIN
MAX_SNP_GAP_BP = m21.MAX_SNP_GAP_BP
GAP_TOL = m21.GAP_TOL
FROH_MIN_MB = m21.FROH_MIN_MB
OUTLIER_F = m21.OUTLIER_F
MIN_KEEP_MB = m21.MIN_KEEP_MB
AF_PRE = m21.AF_PRE
L_GRID = m21.L_GRID
THR_PC = m21.THR_PC
WINDOW_BP = m16.WINDOW_BP

N_BINS = 5
B_BOOT = 500
SEED = 17
OUT = HERE / "per_locus_empirical.tsv"
OUT_TXT = HERE / "per_locus_empirical.txt"


def roh_segments(hom, pos):
    """Like m21.roh_lengths but return (start_bp, end_bp) arrays (gap-tolerant)."""
    m = hom.copy()
    if GAP_TOL > 0:
        pad = np.concatenate(([1], m.astype(np.int8), [1])); dd = np.diff(pad)
        hs = np.flatnonzero(dd == -1); he = np.flatnonzero(dd == 1)
        short = (he - hs) <= GAP_TOL
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
        return np.empty(0, np.int64), np.empty(0, np.int64)
    return pos[starts], pos[ends]


def decode_cum(chrom):
    """Return (rate_per_window, cum_cM, n_windows) for a chrom from the deCODE maps."""
    pat = m16.load_decode_map(m16.PAT_MAP, chrom)
    mat = m16.load_decode_map(m16.MAT_MAP, chrom)
    if not pat:
        return None
    rate, cum = m16.build_cm(pat, mat, max(pat) + WINDOW_BP)
    return rate, cum


def cum_cm_at(cum, rate, bp):
    w = bp // WINDOW_BP
    if w + 1 >= len(cum):
        return cum[-1]
    return cum[w] + rate[w] * (bp - w * WINDOW_BP) / WINDOW_BP


def main():
    t0 = time.time()
    chroms = sys.argv[1:] or [f"chr{n}" for n in range(1, 23)]

    # --- Pass 0: genome-wide window-rate distribution -> bin edges + footprint ---
    all_rates = []
    cuminfo = {}
    for chrom in chroms:
        ci = decode_cum(chrom)
        if ci is None:
            continue
        rate, cum = ci
        cuminfo[chrom] = ci
        all_rates.extend([r for r in rate if r > 0])
    all_rates = np.array(all_rates)
    edges = np.quantile(all_rates, np.linspace(0, 1, N_BINS + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    footprint = np.zeros(N_BINS)          # Mb of autosome per bin
    for chrom in cuminfo:
        rate, _ = cuminfo[chrom]
        for r in rate:
            if r > 0:
                b = int(np.searchsorted(edges, r, side="right") - 1)
                footprint[min(b, N_BINS - 1)] += 1.0   # 1 Mb windows
    bin_mid_rate = [float(np.mean(all_rates[(all_rates >= edges[b]) & (all_rates < edges[b + 1])]))
                    if np.any((all_rates >= edges[b]) & (all_rates < edges[b + 1])) else float("nan")
                    for b in range(N_BINS)]

    # --- Pass 1: parse children, collect per-child (length, bin) segments + F_ROH ---
    sp = m21.load_superpop(); kids = m21.load_children(sp)
    # per child: list of (length_mb, bin_idx); and span, froh-burden
    child_segs = []   # list over children of np arrays: rows (len, bin)
    child_meta = []   # (pop, span_mb, burden_ge1)
    # accumulators keyed by sample name across chroms
    segs_by_kid = {}  # name -> list of (len, bin)
    span_by_kid = {}  # name -> span
    burd_by_kid = {}  # name -> ROH>=1Mb burden
    pop_by_kid = {}

    for chrom in chroms:
        if chrom not in cuminfo:
            continue
        rate, cum = cuminfo[chrom]
        vcf = HERE / "chr22_phased.vcf.gz" if chrom == "chr22" else DATA_DIR_OTHER / f"{chrom}_phased.vcf.gz"
        if not vcf.exists():
            print(f"  [{chrom}] SKIP (no vcf)"); continue
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]; break
        colpop = {p: [] for p in POPS}
        for j, s in enumerate(samples):
            if s in kids:
                colpop[kids[s]].append((j, s))
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
                    rows[p].append(bytes(1 if gts[c][0] == gts[c][2] else 0 for c, _ in cols))
                    pos[p].append(p1)
        for p in POPS:
            if not rows[p]:
                continue
            cols = colpop[p]
            mat = np.frombuffer(b"".join(rows[p]), dtype=np.int8).reshape(len(rows[p]), len(cols))
            pa = np.asarray(pos[p], dtype=np.int64)
            span = (pa[-1] - pa[0]) / 1e6
            for jc, (col, name) in enumerate(cols):
                st, en = roh_segments(mat[:, jc].astype(bool), pa)
                span_by_kid[name] = span_by_kid.get(name, 0.0) + span
                pop_by_kid[name] = p
                lst = segs_by_kid.setdefault(name, [])
                burd = burd_by_kid.get(name, 0.0)
                for s_bp, e_bp in zip(st, en):
                    length = (e_bp - s_bp) / 1e6
                    if length <= MIN_KEEP_MB:
                        continue
                    dcm = cum_cm_at(cum, rate, e_bp) - cum_cm_at(cum, rate, s_bp)
                    mean_r = dcm / length if length > 0 else 0.0
                    b = int(np.searchsorted(edges, mean_r, side="right") - 1)
                    b = min(max(b, 0), N_BINS - 1)
                    lst.append((length, b))
                    if length >= FROH_MIN_MB:
                        burd += length
                burd_by_kid[name] = burd
        print(f"  [{chrom}] parsed ({time.time()-t0:.0f}s)"); sys.stdout.flush()

    # assemble + screen (F_ROH <= OUTLIER_F = the calibration background)
    kept = []
    for name, segs in segs_by_kid.items():
        span = span_by_kid[name]
        froh = burd_by_kid.get(name, 0.0) / span if span else 0.0
        if froh <= OUTLIER_F:
            kept.append(np.array(segs, dtype=float) if segs else np.empty((0, 2)))
    n_kid = len(kept)
    print(f"\n  screened children kept: {n_kid}")

    # --- empirical p_background and L* per bin (+ children bootstrap) ---
    def lstar_for(children_segs, bin_idx):
        segs = np.concatenate([c[c[:, 1] == bin_idx, 0] for c in children_segs if len(c)]) \
            if children_segs else np.empty(0)
        expo = len(children_segs) * footprint[bin_idx]
        if segs.size == 0 or expo == 0:
            return float("inf")
        s = np.sort(segs.astype(np.float64))
        prefix = np.concatenate(([0.0], np.cumsum(s)))
        idx = np.searchsorted(s, L_GRID, side="right")
        emp = (prefix[-1] - prefix[idx] - L_GRID * (s.size - idx)) / expo
        hit = np.flatnonzero(emp <= THR_PC)
        return float(L_GRID[hit[0]]) if hit.size else float("inf")

    rng = np.random.default_rng(SEED)
    point = [lstar_for(kept, b) for b in range(N_BINS)]
    boot = {b: [] for b in range(N_BINS)}
    for _ in range(B_BOOT):
        samp = [kept[i] for i in rng.integers(0, n_kid, n_kid)]
        for b in range(N_BINS):
            boot[b].append(lstar_for(samp, b))
    ci = {b: (float(np.nanpercentile([x for x in boot[b] if np.isfinite(x)] or [np.nan], 2.5)),
              float(np.nanpercentile([x for x in boot[b] if np.isfinite(x)] or [np.nan], 97.5)))
          for b in range(N_BINS)}

    # closed-form L* per bin for contrast (genome-wide mean H-bar, bin mean rate)
    hbar_gw = _genomewide_mean_hbar()
    rows_out = []
    for b in range(N_BINS):
        n_seg = int(sum(int((c[:, 1] == b).sum()) for c in kept if len(c)))
        cf = m16.min_callable_length(bin_mid_rate[b], hbar_gw, m21.PI) if hbar_gw else float("nan")
        rows_out.append((b, edges[b], edges[b + 1], bin_mid_rate[b], footprint[b],
                         n_seg, point[b], ci[b][0], ci[b][1], cf))

    hdr = ("bin\trate_lo\trate_hi\tmean_rate_cMperMb\tfootprint_Mb\tn_segments\t"
           "Lstar_emp_Mb\tCI2.5\tCI97.5\tLstar_closedform_Mb")
    lines = [hdr]
    for r in rows_out:
        lines.append("\t".join(
            (f"{r[0]}", f"{r[1]:.3f}" if np.isfinite(r[1]) else "-inf",
             f"{r[2]:.3f}" if np.isfinite(r[2]) else "inf", f"{r[3]:.3f}", f"{r[4]:.0f}",
             f"{r[5]}", f"{r[6]:.2f}", f"{r[7]:.2f}", f"{r[8]:.2f}", f"{r[9]:.2f}")))
    OUT.write_text("\n".join(lines) + "\n")

    txt = [f"Per-locus EMPIRICAL decisive length (screened children={n_kid}, "
           f"pooled; {N_BINS} recombination-rate bins, B={B_BOOT} bootstrap)",
           f"chroms={','.join(chroms)}  THR_PC={THR_PC:.2e}  prior=1st-cousin",
           "=" * 74,
           f"{'rate cM/Mb (bin mean)':>22} {'n_seg':>7} {'L* empirical (95% CI)':>26} {'L* closed-form':>15}"]
    for r in rows_out:
        lo = f"{r[7]:.2f}"; hi = f"{r[8]:.2f}"
        emp = "inf" if not np.isfinite(r[6]) else f"{r[6]:.2f}"
        txt.append(f"{r[3]:>22.2f} {r[5]:>7d} {emp:>10} ({lo}-{hi}) Mb     {r[9]:>10.2f} Mb")
    direction = ("DECREASES with rate (thesis supported)"
                 if np.isfinite(point[0]) and np.isfinite(point[-1]) and point[-1] < point[0]
                 else "non-monotone / inconclusive (check tail counts)")
    txt += ["", f"Empirical L*: {direction}.",
            "Empirical vs closed-form: the closed form is the analytic prediction; "
            "divergence quantifies its anti-conservatism per bin."]
    OUT_TXT.write_text("\n".join(txt) + "\n")
    print("\n".join(txt))
    print(f"\nwrote {OUT.name}, {OUT_TXT.name}  ({time.time()-t0:.0f}s)")


def _genomewide_mean_hbar():
    div = HERE / "cross_pop_hap_diversity.tsv"
    if not div.exists():
        return None
    vals = []
    with div.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals.append(float(f[ix["Hbar"]]))
    return float(np.mean(vals)) if vals else None


if __name__ == "__main__":
    main()
