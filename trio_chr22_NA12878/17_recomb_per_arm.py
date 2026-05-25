"""
Recombination summary per chromosome ARM (p / q), from the deCODE / Palsson 2024
sex-specific crossover maps (GRCh38; cMperMb at 1 Mb window centers).

This is a pure map summary -- it reads no genotypes and is independent of the
ROH / H-bar pipeline. The cMperMb column is genuinely sex-specific between the
pat and mat files (e.g. chr1:500 kb is paternal 1.32 vs maternal 0.05), so the
paternal vs maternal arm totals below are real, not duplicated.

Per autosomal arm it reports:
  - covered physical span (Mb) and number of 1 Mb map windows
  - genetic length in cM: paternal, maternal, sex-averaged  (sum of cMperMb,
    since each window is 1 Mb wide so cM == sum of per-Mb rates)
  - mean recombination rate cM/Mb (sex-averaged)
  - Haldane recombination fraction theta across the whole arm (sex-averaged):
        theta = 0.5 * (1 - exp(-2d)),   d = cM / 100  Morgans   [no interference]
    theta saturates toward 0.5 for long arms -- it is the chance two loci at the
    arm's ends end up in different gametes.

Arm split uses GRCh38 centromere spans (UCSC cytoBand 'acen', merged). Windows
whose center falls inside the centromere span are pericentromeric and excluded
from both arms (recombination there is strongly suppressed, so this barely moves
arm totals). Acrocentric p-arms (13p,14p,15p,21p,22p) have no map coverage and
are reported as NA.

Inputs (already on disk): external/.../maps.pat.tsv, maps.mat.tsv
Output: recomb_per_arm.tsv   (+ printed summary incl. genome-wide p/q and M:F)
Usage:  python 17_recomb_per_arm.py
"""

import math
from pathlib import Path

HERE = Path(__file__).parent
DECODE_DIR = HERE / "external" / "palsson2024_deCODE_maps" / "DecodeGenetics-PalssonEtAl_Nature_2024-8e49794" / "data" / "maps"
PAT_MAP = DECODE_DIR / "maps.pat.tsv"
MAT_MAP = DECODE_DIR / "maps.mat.tsv"
OUT = HERE / "recomb_per_arm.tsv"

WINDOW_BP = 1_000_000

# GRCh38 centromere spans (bp), UCSC cytoBand 'acen' merged. Only ~Mb precision
# is needed: pericentromeric windows carry little cM.
CENTROMERE = {
    "chr1": (122026459, 125184587),  "chr2": (92188145, 94090557),
    "chr3": (90772458, 93655574),    "chr4": (49712061, 51743951),
    "chr5": (46485900, 50059807),    "chr6": (58553888, 59829934),
    "chr7": (58169653, 61528020),    "chr8": (44033744, 45877265),
    "chr9": (43389635, 45518558),    "chr10": (39686682, 41593521),
    "chr11": (51078348, 54425074),   "chr12": (34769407, 37185252),
    "chr13": (16000000, 18051248),   "chr14": (16000000, 18173523),
    "chr15": (17083673, 19725254),   "chr16": (36311158, 38280682),
    "chr17": (22813679, 26885980),   "chr18": (15460899, 20861206),
    "chr19": (24498980, 27190874),   "chr20": (26436232, 30038348),
    "chr21": (10864560, 12915808),   "chr22": (12954788, 15054318),
}
AUTOSOMES = [f"chr{n}" for n in range(1, 23)]


def load_cmpermb(path):
    """{chrom: {pos_center: cMperMb}} from a deCODE map file (col0 Chr, col1 pos,
    col3 cMperMb)."""
    out = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Chr"):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                out.setdefault(f[0], {})[int(f[1])] = float(f[3])
            except (IndexError, ValueError):
                continue
    return out


def haldane_theta(cm):
    return 0.5 * (1.0 - math.exp(-2.0 * cm / 100.0))


def main():
    pat = load_cmpermb(PAT_MAP)
    mat = load_cmpermb(MAT_MAP)

    rows = []
    tot = {"p": [0.0, 0.0], "q": [0.0, 0.0]}  # arm -> [pat_cM, mat_cM]
    for chrom in AUTOSOMES:
        cstart, cend = CENTROMERE[chrom]
        arms = {a: {"n": 0, "pat": 0.0, "mat": 0.0, "lo": None, "hi": None}
                for a in ("p", "q")}
        for pos in sorted(pat.get(chrom, {})):
            if cstart <= pos <= cend:
                continue  # pericentromeric window -> neither arm
            a = arms["p" if pos < cstart else "q"]
            a["n"] += 1
            a["pat"] += pat[chrom][pos]
            a["mat"] += mat.get(chrom, {}).get(pos, pat[chrom][pos])
            a["lo"] = pos if a["lo"] is None else min(a["lo"], pos)
            a["hi"] = pos if a["hi"] is None else max(a["hi"], pos)

        for arm in ("p", "q"):
            a = arms[arm]
            if a["n"] == 0:
                rows.append((chrom, arm, "NA", 0, "NA", "NA", "NA", "NA", "NA"))
                continue
            cm_avg = 0.5 * (a["pat"] + a["mat"])
            span_mb = (a["hi"] - a["lo"]) / 1e6 + 1.0  # window-center span + 1 Mb
            tot[arm][0] += a["pat"]
            tot[arm][1] += a["mat"]
            rows.append((chrom, arm, f"{span_mb:.1f}", a["n"],
                         f"{cm_avg:.2f}", f"{a['pat']:.2f}", f"{a['mat']:.2f}",
                         f"{cm_avg / a['n']:.3f}", f"{haldane_theta(cm_avg):.4f}"))

    hdr = ("chrom", "arm", "span_Mb", "n_win", "cM_sexavg", "cM_pat", "cM_mat",
           "cMperMb_sexavg", "theta_haldane_sexavg")
    with OUT.open("w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"  deCODE / Palsson 2024 GRCh38 maps -> {OUT}\n")
    print("  " + "\t".join(hdr))
    for r in rows:
        print("  " + "\t".join(str(x) for x in r))

    p_avg = 0.5 * (tot["p"][0] + tot["p"][1])
    q_avg = 0.5 * (tot["q"][0] + tot["q"][1])
    pat_tot = tot["p"][0] + tot["q"][0]
    mat_tot = tot["p"][1] + tot["q"][1]
    print("\n  genome-wide autosomal genetic length (sex-averaged):")
    print(f"    p-arms {p_avg:.1f} cM | q-arms {q_avg:.1f} cM | total {p_avg + q_avg:.1f} cM")
    print(f"  paternal total {pat_tot:.1f} cM | maternal total {mat_tot:.1f} cM | "
          f"male:female ratio {pat_tot / mat_tot:.3f} (expect <1; females recombine more)")


if __name__ == "__main__":
    main()
