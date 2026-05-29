"""
Cleaner non-IBD null from the 1000G TRIO CHILDREN (addresses validation leakage /
artificial-pair / cryptic-relatedness criticisms).

The 602 trio children are the "related" set, DISJOINT from the 2,504 unrelated
panel used to estimate H-bar. A nominally-outbred child's two homologs ARE the
real chance/background of an outbred individual -> the non-IBD null the posterior
must beat, measured in real people not used to build the model.

This version adds a per-child cryptic-relatedness / endogamy screen: each child's
genome-wide ROH burden (sum of ROH >= 1 Mb / scanned length = F_ROH) is computed,
and children with F_ROH above OUTLIER_F (recent shared ancestry) are flagged. The
background null is reported BOTH including all children and excluding flagged ones,
so true population background is separated from a few related samples.

Prior-free output: Bayes factor BF(L) = c / p_background(L) (weight of evidence).

Outputs: trio_null_summary.txt, trio_null_pchance.tsv, trio_null_perchild.tsv
Usage:  python 21_trio_background_null.py [chrom ...]   (no args = all 22 autosomes)
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
OUT_SUM = HERE / "trio_null_summary.txt"
OUT_PC = HERE / "trio_null_pchance.tsv"
OUT_CHILD = HERE / "trio_null_perchild.tsv"

POPS = ["EUR", "AFR", "EAS", "SAS", "AMR"]
MAF_MIN = 0.05
GAP_TOL = 1
MAX_SNP_GAP_BP = 1_000_000
MIN_KEEP_MB = 0.005
FROH_MIN_MB = 1.0                 # ROH >= this count toward the F_ROH burden
OUTLIER_F = 0.0156                # F_ROH above this (~2nd-cousin) flags recent shared ancestry
WINDOW_BP = 1_000_000
BLOCK_CM = 0.5
H_FLOOR = 1e-4
GENO_ERR = 0.001
PI = 0.0625
T_DEC = 0.95
C_IBD = (1 - GENO_ERR) ** 1000
L_GRID = np.round(np.arange(0.1, 12.001, 0.05), 3)
L_REPORT = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
AF_PRE = {p: f"AF_{p}=" for p in POPS}
THR_PC = PI * C_IBD * (1 - T_DEC) / (T_DEC * (1 - PI))


def load_superpop():
    sp = {}
    with PANEL.open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                sp[f[0]] = f[2]
    return sp


def load_children(sp):
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
            p = f[ix["population"]]
            if p in POPS:
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


def parse_children(chroms):
    """Parse trio children per population -> (childsegs, childtot, childspan,
    nchild). childsegs[p][j] = list of per-chrom ROH-segment arrays for child j;
    childtot[p][j] = ROH>=FROH_MIN_MB burden (Mb); childspan[p][j] = scanned span
    (Mb). Extracted from main() so the bootstrap (script 24) reuses one parse."""
    t0 = time.time()
    sp = load_superpop()
    kids = load_children(sp)
    childsegs = {p: None for p in POPS}   # per child: list of seg arrays
    childtot = {p: None for p in POPS}    # per child: ROH>=FROH_MIN_MB burden (Mb)
    childspan = {p: None for p in POPS}   # per child: scanned span (Mb)
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
                sl = roh_lengths(mat[:, j].astype(bool), pa)
                childtot[p][j] += float(sl[sl >= FROH_MIN_MB].sum())
                childspan[p][j] += span
                slk = sl[sl > MIN_KEEP_MB]
                if slk.size:
                    childsegs[p][j].append(slk)
        print(f"  [{chrom}] children/pop " + " ".join(f"{p}:{nchild[p]}" for p in POPS)
              + f"  ({time.time()-t0:.0f}s)")
        sys.stdout.flush()
    return childsegs, childtot, childspan, nchild


def main():
    t0 = time.time()
    chroms = sys.argv[1:] or [f"chr{n}" for n in range(1, 23)]
    childsegs, childtot, childspan, nchild = parse_children(chroms)

    div = load_div(chroms[0])
    rmean = np.nanmean([v[1] for p in POPS for v in div[p].values()]) if any(div.values()) else 1.2
    hbar_gw = {p: (np.nanmean([v[0] for v in div[p].values()]) if div[p] else float("nan")) for p in POPS}

    # per-child F_ROH + outlier flags; both nulls
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
        emp_all[p], Lall[p] = agg_null(childsegs[p], childspan[p])
        keep = [j for j in range(nchild[p]) if not outl[j]]
        emp_clean[p], Lclean[p] = agg_null([childsegs[p][j] for j in keep],
                                           [childspan[p][j] for j in keep])

    with OUT_CHILD.open("w", encoding="utf-8") as fh:
        fh.write("pop\tchild_idx\tROH_ge1Mb_Mb\tspan_Mb\tF_ROH\tflag\n")
        for r in per_child_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with OUT_PC.open("w", encoding="utf-8") as fh:
        cols = [p for p in POPS if p in emp_clean and emp_clean[p] is not None]
        fh.write("L_Mb\t" + "\t".join(f"{p}_pchance_clean\t{p}_log10BF_clean" for p in cols) + "\n")
        for k, L in enumerate(L_GRID):
            cell = []
            for p in cols:
                pc = max(emp_clean[p][k], 1e-12); cell.append(f"{emp_clean[p][k]:.3e}\t{math.log10(C_IBD/pc):.2f}")
            fh.write(f"{L:.2f}\t" + "\t".join(cell) + "\n")

    with OUT_SUM.open("w", encoding="utf-8") as fh:
        fh.write("# Trio-children background null (leakage-free, cryptic-relatedness-screened)\n")
        fh.write(f"# chroms={','.join(chroms)} children={sum(nchild.values())} "
                 f"PI={PI} GAP_TOL={GAP_TOL} OUTLIER_F={OUTLIER_F} (ROH>={FROH_MIN_MB}Mb)\n")
        fh.write(f"# prior-free Bayes factor BF(L)=c/p_background(L), c={C_IBD:.3f}. "
                 f"wall={time.time()-t0:.0f}s\n\n")
        fh.write("Minimum 'callable' length L* (posterior>={:.2f}, pi={}) per population:\n".format(T_DEC, PI))
        fh.write("pop\tn_children\tn_outliers(F_ROH>%.4f)\tL*_all\tL*_screened\tHbar_gw\tanalytic_L*\n" % OUTLIER_F)
        for p in POPS:
            if p not in Lall:
                continue
            b = max(hbar_gw[p], H_FLOOR)
            anaL = (max(math.log(THR_PC) / math.log(b) * BLOCK_CM / rmean, BLOCK_CM / rmean)
                    if b < 1 and rmean > 0 else float("nan"))
            fh.write(f"{p}\t{nchild[p]}\t{nout[p]}\t{Lall[p]:.2f}\t{Lclean[p]:.2f}\t"
                     f"{hbar_gw[p]:.4f}\t{anaL:.2f}\n")
        fh.write("\nWeight of evidence log10 BF(L) at tabulated lengths "
                 "(cryptic-relatedness-screened background):\n")
        cols = [p for p in POPS if p in emp_clean and emp_clean[p] is not None]
        fh.write("L_Mb\t" + "\t".join(cols) + "\n")
        for L in L_REPORT:
            k = int(np.argmin(np.abs(L_GRID - L)))
            fh.write(f"{L}\t" + "\t".join(
                f"{math.log10(C_IBD/max(emp_clean[p][k],1e-12)):.2f}" for p in cols) + "\n")
        fh.write(f"\nthr p_chance for posterior>={T_DEC}: {THR_PC:.2e}. Children are DISJOINT "
                 "from the unrelated panel used to estimate H-bar. L*_screened excludes "
                 "children with recent shared ancestry (F_ROH outliers).\n")

    print(f"\n  -> {OUT_SUM}\n  -> {OUT_PC}\n  -> {OUT_CHILD}\n  total {time.time()-t0:.0f}s")
    for p in POPS:
        if p in Lall:
            print(f"    {p}: L*_all={Lall[p]:.2f}  L*_screened={Lclean[p]:.2f} Mb  "
                  f"(children={nchild[p]}, outliers={nout[p]})")


if __name__ == "__main__":
    main()
