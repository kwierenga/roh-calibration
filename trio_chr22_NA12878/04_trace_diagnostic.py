"""
Diagnostic: bin the chr22 trace into 250 kb windows along the chromosome.
For each window and each parent, count how many informative sites called
haplotype 0 vs haplotype 1. Output a TSV and an ASCII chromosome trace
so we can SEE where the haplotype state actually changes.

This lets us distinguish:
  - Real crossover: a clean transition (window state goes 0,0,0,0 -> 1,1,1,1)
  - Phasing noise: ragged windows with mixed 0/1 counts
  - No crossover: uniform state across the chromosome
"""

import gzip
from pathlib import Path

HERE = Path(__file__).parent
IN_TSV = HERE / "trio_chr22.tsv.gz"
OUT_TSV = HERE / "trace_chr22_250kb.tsv"
OUT_ASCII = HERE / "trace_chr22_ascii.txt"

WINDOW_BP = 250_000


def parse_gt(gt):
    if "|" not in gt:
        return None
    a, b = gt.split("|", 1)
    try:
        return int(a), int(b)
    except ValueError:
        return None


def main():
    # bin counts: windows[start_window_index] = [pat0, pat1, mat0, mat1]
    windows = {}

    with gzip.open(IN_TSV, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pos = int(f[idx["POS"]])
            father = parse_gt(f[idx["FATHER"]])
            mother = parse_gt(f[idx["MOTHER"]])
            child = parse_gt(f[idx["CHILD"]])
            if father is None or mother is None or child is None:
                continue

            wstart = (pos // WINDOW_BP) * WINDOW_BP
            w = windows.setdefault(wstart, [0, 0, 0, 0])

            # paternal: only informative if father is het
            if father[0] != father[1]:
                # child[0] is paternally-inherited under convention A
                ca = child[0]
                if ca == father[0]:
                    w[0] += 1  # pat hap 0
                elif ca == father[1]:
                    w[1] += 1  # pat hap 1

            # maternal: only informative if mother is het
            if mother[0] != mother[1]:
                ca = child[1]
                if ca == mother[0]:
                    w[2] += 1  # mat hap 0
                elif ca == mother[1]:
                    w[3] += 1  # mat hap 1

    # write TSV
    with OUT_TSV.open("w") as out:
        out.write("window_start\twindow_end\tpat_hap0\tpat_hap1\tpat_dominant\tpat_purity\tmat_hap0\tmat_hap1\tmat_dominant\tmat_purity\n")
        for wstart in sorted(windows):
            p0, p1, m0, m1 = windows[wstart]
            wend = wstart + WINDOW_BP
            ptot = p0 + p1
            mtot = m0 + m1
            pdom = "0" if p0 > p1 else ("1" if p1 > p0 else "-")
            mdom = "0" if m0 > m1 else ("1" if m1 > m0 else "-")
            ppur = max(p0, p1) / ptot if ptot else 0
            mpur = max(m0, m1) / mtot if mtot else 0
            out.write(f"{wstart}\t{wend}\t{p0}\t{p1}\t{pdom}\t{ppur:.3f}\t{m0}\t{m1}\t{mdom}\t{mpur:.3f}\n")

    # ASCII chromosome trace
    # Each window becomes one character. We use:
    #   '0' if hap0 dominates with purity >= 0.9
    #   '1' if hap1 dominates with purity >= 0.9
    #   'x' if mixed (purity 0.5-0.9)
    #   '.' if no informative sites
    def code(hap0, hap1):
        tot = hap0 + hap1
        if tot < 5:
            return "."
        pur = max(hap0, hap1) / tot
        if pur >= 0.9:
            return "0" if hap0 > hap1 else "1"
        elif pur >= 0.7:
            return "0" if hap0 > hap1 else "1"  # still dominant, but flag with lowercase below
        else:
            return "x"

    def code_with_purity_flag(hap0, hap1):
        tot = hap0 + hap1
        if tot < 5:
            return "."
        pur = max(hap0, hap1) / tot
        dom = "0" if hap0 > hap1 else "1"
        if pur >= 0.95:
            return dom            # very pure
        elif pur >= 0.80:
            return dom.lower() if dom.isalpha() else dom  # somewhat pure
        else:
            return "x"            # mixed

    sorted_starts = sorted(windows)
    pat_line = "".join(code_with_purity_flag(windows[w][0], windows[w][1]) for w in sorted_starts)
    mat_line = "".join(code_with_purity_flag(windows[w][2], windows[w][3]) for w in sorted_starts)

    # also a coordinate ruler every 10 Mb
    ruler = []
    for w in sorted_starts:
        mb = w // 1_000_000
        if w % 10_000_000 == 0:
            ruler.append(str(mb // 10))
        elif w % 1_000_000 == 0:
            ruler.append(str(mb % 10))
        else:
            ruler.append(" ")
    ruler_line = "".join(ruler)

    with OUT_ASCII.open("w") as out:
        out.write(f"# chr22 haplotype-state trace in {WINDOW_BP/1000:.0f} kb windows\n")
        out.write("# '0' or '1' = haplotype state, dominant > 95% pure\n")
        out.write("# 'x' = mixed (< 80% pure)\n")
        out.write("# '.' = fewer than 5 informative sites\n")
        out.write(f"# total windows: {len(sorted_starts)}\n")
        out.write(f"# first window start: {sorted_starts[0]:,}\n")
        out.write(f"# last window start:  {sorted_starts[-1]:,}\n\n")

        # wrap at 100 chars per line for readability
        WRAP = 100
        n = len(pat_line)
        for chunk_start in range(0, n, WRAP):
            chunk_end = min(chunk_start + WRAP, n)
            pos_start = sorted_starts[chunk_start]
            pos_end = sorted_starts[chunk_end - 1] + WINDOW_BP
            out.write(f"\nchr22:{pos_start:,}-{pos_end:,}\n")
            out.write(f"  ruler: {ruler_line[chunk_start:chunk_end]}\n")
            out.write(f"  pat:   {pat_line[chunk_start:chunk_end]}\n")
            out.write(f"  mat:   {mat_line[chunk_start:chunk_end]}\n")

    print(f"  windows analyzed: {len(sorted_starts)}")
    print(f"  -> {OUT_TSV}")
    print(f"  -> {OUT_ASCII}")
    print()
    print("CHROMOSOME-WIDE TRACE (one char per 250 kb window):")
    print("  '0'/'1' = haplotype state (>= 95% pure); 'x' = mixed; '.' = sparse")
    print()
    print(f"  pat: {pat_line}")
    print(f"  mat: {mat_line}")


if __name__ == "__main__":
    main()
