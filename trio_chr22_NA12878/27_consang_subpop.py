"""
Per-SAS-SUBPOP background null (PJL/GIH/STU/BEB/ITU) — addresses reviewer #2's
"the screen captures population background, not just cryptic relatedness" concern
by breaking out 1000G SAS into its 5 source panels, of which PJL (Punjabi in Lahore,
N=96) is the highest-consanguinity (~30-40% in source) and most clinically relevant.

Identical methodology to 21_trio_background_null.py but operating at SAS-subpop
level. Uses AF_SAS (superpop-level) as the MAF filter for all 5 subpops since
1000G HC VCFs do not carry subpop-level AF fields. H-bar is taken at SAS superpop
level from cross_pop_hap_diversity.tsv (slight approximation; per-subpop H-bar
would refine this further but is orthogonal to the L*_all vs L*_screened question).

Output: trio_null_consang_summary.txt, trio_null_consang_pchance.tsv,
trio_null_consang_perchild.tsv
Usage: python 27_consang_subpop.py [chrom ...]   (no args = all 22 autosomes)
"""

import gzip
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DATA_DIR_OTHER = HERE / "all_autosomes"
PED = HERE / "pedigree.txt"
PANEL = HERE / "samples_2504_pop.panel"
DIVTSV = HERE / "cross_pop_hap_diversity.tsv"
OUT_SUM = HERE / "trio_null_consang_summary.txt"
OUT_PC = HERE / "trio_null_consang_pchance.tsv"
OUT_CHILD = HERE / "trio_null_consang_perchild.tsv"

POPS = ["PJL", "GIH", "STU", "BEB", "ITU"]
SUPERPOP_OF = {p: "SAS" for p in POPS}
MAF_MIN = 0.05
GAP_TOL = 1
MAX_SNP_GAP_BP = 1_000_000
MIN_KEEP_MB = 0.005
FROH_MIN_MB = 1.0
OUTLIER_F = 0.0156
WINDOW_BP = 1_000_000
BLOCK_CM = 0.5
H_FLOOR = 1e-4
GENO_ERR = 0.001
PI = 0.0625
T_DEC = 0.95
C_IBD = (1 - GENO_ERR) ** 1000
L_GRID = np.round(np.arange(0.1, 12.001, 0.05), 3)
L_REPORT = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
AF_PRE = "AF_SAS="
THR_PC = PI * C_IBD * (1 - T_DEC) / (T_DEC * (1 - PI))


def load_subpop():
    """sample -> SAS subpop (only for samples in the panel that are SAS)."""
    sp = {}
    with PANEL.open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] == "SAS":
                sp[f[0]] = f[1]
    return sp


def load_children(sp):
    """Inherit subpop from either parent (within SAS only)."""
    kids = {}
    with PED.open() as fh:
        next(fh)
        for line in fh:
            f = line.split()
            if len(f) >= 3 and f[1] != "0" and f[2] != "0":
                pop = sp.get(f[1]) or sp.get(f[2])
                if pop in POPS:
                    kids[f[0]] = pop
    return kids


def load_div(chrom):
    out = {p: {} for p in POPS}
    if not DIVTSV.exists():
        return out
    with DIVTSV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["chrom"]] != chrom:
                continue
            sup = f[ix["population"]]
            if sup == "SAS":
                for p in POPS:
                    out[p][int(f[ix["window_start"]])] = (float(f[ix["Hbar"]]), float(f[ix["cMperMb"]]))
    return out


def roh_lengths(hom, pos):
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
        return np.empty(0, np.float32)
    return ((pos[ends] - pos[starts]) / 1e6).astype(np.float32)


def agg_null(seg_lists, spans):
    segs = [a for cs in seg_lists for a in cs]
    expo = float(sum(spans))
    if not segs or expo == 0:
        return None, float("nan")
    s = np.sort(np.concatenate(segs).astype(np.float64))
    prefix = np.concatenate(([0.0], np.cumsum(s)))
    idx = np.searchsorted(s, L_GRID, side="right")
    emp = (prefix[-1] - prefix[idx] - L_GRID * (s.size - idx)) / expo
    hit = np.flatnonzero(emp <= THR_PC)
    return emp, (float(L_GRID[hit[0]]) if hit.size else float("inf"))


def main():
    t0 = time.time()
    chroms = sys.argv[1:] or [f"chr{n}" for n in range(1, 23)]
    sp = load_subpop()
    kids = load_children(sp)
    childsegs = {p: None for p in POPS}
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
            if childsegs[p] is None and colpop[p]:
                n = len(colpop[p])
                childsegs[p] = [[] for _ in range(n)]
                childtot[p] = [0.0] * n; childspan[p] = [0.0] * n; nchild[p] = n
        rows_buf = []; pos_buf = []
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
                    continue
                # SAS-superpop MAF filter applied uniformly to all 5 subpops
                af_sas = None
                for kv in f[7].split(";"):
                    if kv.startswith(AF_PRE):
                        try:
                            af_sas = float(kv[len(AF_PRE):])
                        except ValueError:
                            af_sas = None
                        break
                if af_sas is None or min(af_sas, 1 - af_sas) < MAF_MIN:
                    continue
                gts = f[9:]; p1 = int(f[1])
                # build hom row across ALL subpop columns (union)
                all_cols = []
                for p in POPS:
                    all_cols.extend(colpop[p])
                if not all_cols:
                    continue
                rows_buf.append(bytes(1 if gts[c][0] == gts[c][2] else 0 for c in all_cols))
                pos_buf.append(p1)
        if not rows_buf:
            continue
        # split the union-row matrix back per subpop
        all_cols_idx = []; subpop_slices = {}
        cur = 0
        for p in POPS:
            n_p = len(colpop[p])
            subpop_slices[p] = (cur, cur + n_p)
            cur += n_p
            all_cols_idx.extend(colpop[p])
        ncol_union = cur
        mat_union = np.frombuffer(b"".join(rows_buf), dtype=np.int8).reshape(len(rows_buf), ncol_union)
        pa = np.asarray(pos_buf, dtype=np.int64)
        span = (pa[-1] - pa[0]) / 1e6
        for p in POPS:
            s_lo, s_hi = subpop_slices[p]
            if s_hi == s_lo:
                continue
            for j in range(s_lo, s_hi):
                sl = roh_lengths(mat_union[:, j].astype(bool), pa)
                local_j = j - s_lo
                childtot[p][local_j] += float(sl[sl >= FROH_MIN_MB].sum())
                childspan[p][local_j] += span
                slk = sl[sl > MIN_KEEP_MB]
                if slk.size:
                    childsegs[p][local_j].append(slk)
        print(f"  [{chrom}] children/subpop " + " ".join(f"{p}:{nchild[p]}" for p in POPS)
              + f"  ({time.time()-t0:.0f}s)")
        sys.stdout.flush()

    div = load_div(chroms[0])
    rmean = np.nanmean([v[1] for p in POPS for v in div[p].values()]) if any(div.values()) else 1.2
    hbar_gw = {p: (np.nanmean([v[0] for v in div[p].values()]) if div[p] else float("nan")) for p in POPS}

    per_child_rows = []
    emp_all = {}; emp_clean = {}; Lall = {}; Lclean = {}; nout = {}
    for p in POPS:
        if childsegs[p] is None:
            continue
        froh = [childtot[p][j] / childspan[p][j] if childspan[p][j] else 0.0
                for j in range(nchild[p])]
        outl = [froh[j] > OUTLIER_F for j in range(nchild[p])]
        nout[p] = sum(outl)
        for j in range(nchild[p]):
            per_child_rows.append((p, j, round(childtot[p][j], 2), round(childspan[p][j], 1),
                                   round(froh[j], 5), "outlier" if outl[j] else "ok"))
        keep_all = list(range(nchild[p]))
        keep_clean = [j for j in range(nchild[p]) if not outl[j]]
        emp_all[p], Lall[p] = agg_null([childsegs[p][j] for j in keep_all],
                                       [childspan[p][j] for j in keep_all])
        emp_clean[p], Lclean[p] = agg_null([childsegs[p][j] for j in keep_clean],
                                           [childspan[p][j] for j in keep_clean])

    with OUT_CHILD.open("w", encoding="utf-8") as fh:
        fh.write("subpop\tchild_idx\tROH_ge1Mb_Mb\tspan_Mb\tF_ROH\tflag\n")
        for r in per_child_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with OUT_PC.open("w", encoding="utf-8") as fh:
        cols = [p for p in POPS if p in emp_clean and emp_clean[p] is not None]
        fh.write("L_Mb\t"
                 + "\t".join(f"{p}_pchance_all\t{p}_pchance_clean\t{p}_log10BF_clean" for p in cols)
                 + "\n")
        for k, L in enumerate(L_GRID):
            cell = []
            for p in cols:
                pall = emp_all[p][k] if emp_all[p] is not None else float("nan")
                pcle = emp_clean[p][k] if emp_clean[p] is not None else float("nan")
                bf = math.log10(C_IBD / max(pcle, 1e-12)) if pcle == pcle else float("nan")
                cell.append(f"{pall:.6e}\t{pcle:.6e}\t{bf:.3f}")
            fh.write(f"{L}\t" + "\t".join(cell) + "\n")

    with OUT_SUM.open("w", encoding="utf-8") as fh:
        fh.write("# Per-SAS-subpop background null (PJL/GIH/STU/BEB/ITU)\n")
        fh.write(f"# chroms={','.join(chroms)} children={sum(nchild.values())} "
                 f"PI={PI} GAP_TOL={GAP_TOL} OUTLIER_F={OUTLIER_F} (ROH>={FROH_MIN_MB}Mb)\n")
        fh.write(f"# prior-free Bayes factor BF(L)=c/p_background(L), c={C_IBD:.3f}. "
                 f"wall={time.time()-t0:.0f}s\n\n")
        fh.write("Minimum 'callable' length L* (posterior>={:.2f}, pi={}) per SAS subpop:\n".format(T_DEC, PI))
        fh.write("subpop\tn_children\tn_outliers(F_ROH>%.4f)\toutlier_rate\tL*_all\tL*_screened\tdelta_Mb\tHbar_gw_SAS\n"
                 % OUTLIER_F)
        for p in POPS:
            if p not in Lall:
                continue
            rate = nout[p] / nchild[p] if nchild[p] else 0.0
            delta = (Lall[p] - Lclean[p]) if (Lall[p] != float("inf") and Lclean[p] != float("inf")) else float("inf")
            fh.write(f"{p}\t{nchild[p]}\t{nout[p]}\t{rate:.2%}\t{Lall[p]:.2f}\t{Lclean[p]:.2f}\t{delta:.2f}\t"
                     f"{hbar_gw[p]:.4f}\n")
        fh.write("\nWeight of evidence log10 BF(L) at tabulated lengths (screened background):\n")
        cols = [p for p in POPS if p in emp_clean and emp_clean[p] is not None]
        fh.write("L_Mb\t" + "\t".join(cols) + "\n")
        for L in L_REPORT:
            k = int(np.argmin(np.abs(L_GRID - L)))
            fh.write(f"{L}\t" + "\t".join(
                f"{math.log10(C_IBD / max(emp_clean[p][k], 1e-12)):.2f}" for p in cols) + "\n")
        fh.write("\nUNSCREENED log10 BF(L) at tabulated lengths (no F_ROH outlier removal):\n")
        fh.write("L_Mb\t" + "\t".join(cols) + "\n")
        for L in L_REPORT:
            k = int(np.argmin(np.abs(L_GRID - L)))
            fh.write(f"{L}\t" + "\t".join(
                f"{math.log10(C_IBD / max(emp_all[p][k], 1e-12)):.2f}" for p in cols) + "\n")
        fh.write(f"\nthr p_chance for posterior>={T_DEC}: {THR_PC:.2e}. "
                 "MAF filter applied uniformly via AF_SAS for all 5 subpops; "
                 "H-bar approximated at SAS superpop level. "
                 "PJL (Punjabi in Lahore) is the highest-consanguinity 1000G panel.\n")

    print(f"\n  -> {OUT_SUM}\n  -> {OUT_PC}\n  -> {OUT_CHILD}\n  total {time.time()-t0:.0f}s")
    for p in POPS:
        if p in Lall:
            print(f"    {p}: L*_all={Lall[p]:.2f}  L*_screened={Lclean[p]:.2f} Mb  "
                  f"(children={nchild[p]}, outliers={nout[p]}, rate={nout[p]/nchild[p]:.0%})")


if __name__ == "__main__":
    main()
