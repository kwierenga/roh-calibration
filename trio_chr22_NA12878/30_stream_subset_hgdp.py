"""
30_stream_subset_hgdp.py - stream gnomAD v3.1.2 HGDP+1KG VCFs from the public
GCS bucket, filter on the fly, and write a small per-chrom subset to local
disk. The raw .bgz never lands on disk.

Why: the laptop has ~234 GB free; raw all-autosomes is ~1.5-6 TB. Standard
"download-then-filter" doesn't fit. This pipeline keeps disk peak at ~50-200
MB per chrom (the subset, what we keep), ~1-4 GB total across 22 autosomes.

Subset content (matches script 29 output produced from local .bgz):
  - samples : HGDP only (~828 of ~4,150 in the joint callset)
  - variants: biallelic SNVs, single-base REF and ALT, no commas in ALT
  - filter  : MAF (in HGDP subset) >= MAF_MIN
  - columns : INFO stripped to AF_HGDP only; FORMAT stripped to GT

Resumable: skips chroms whose final .vcf.gz already exists. Atomic via
.tmp -> rename so an interrupted chrom never leaves a half-finished output.
Network failures are retried per chrom with exponential backoff.

Usage:
  python 30_stream_subset_hgdp.py            # all 22 autosomes, smallest first
  python 30_stream_subset_hgdp.py 22 21 20   # explicit chroms, given order

Notes:
  - bgzipped raw stays in the gnomAD bucket - we never have a local copy.
  - If chr22.hgdp_tgp.vcf.bgz is already on disk (56 GB), delete it after
    this script produces chr22.hgdp_subset.vcf.gz to reclaim the space.
  - Bandwidth on the wire is still ~the raw size (no way to filter at the
    bucket); the win is purely on disk and on safety (no half-TB temp files).
"""

import gzip
import socket
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
        "release/3.1.2/vcf/genomes/")
FILE = "gnomad.genomes.v3.1.2.hgdp_tgp.chr{c}.vcf.bgz"

HERE = Path(__file__).parent
DDIR = HERE / "hgdp_tgp"
DDIR.mkdir(exist_ok=True)

MAF_MIN = 0.05
N_RETRIES = 3
RETRY_BACKOFF_S = 60
PROGRESS_EVERY_S = 30

# Smallest-first chrom order at gnomAD v3.1.2 (by raw .bgz size). Fast wins
# first, chr1/chr2 last when the pipeline is proven.
SMALLEST_FIRST = ["22", "21", "20", "19", "18", "17", "16", "15", "14",
                  "13", "12", "11", "10", "9", "8", "7", "6", "5", "4",
                  "3", "2", "1"]


def _filter_stream(gz, fout, chrom):
    """Iterate decompressed VCF lines, filter, write minimal subset."""
    n_in = n_pass = n_kept = 0
    hgdp_cols: list[int] = []
    last_print = time.time()
    bytes_seen = 0

    for raw in gz:
        bytes_seen += len(raw)
        line = raw.decode("ascii", "replace")

        if line.startswith("##"):
            if (line.startswith("##fileformat")
                    or line.startswith(f"##contig=<ID=chr{chrom},")
                    or line.startswith("##reference")):
                fout.write(line)
            continue

        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            samples = cols[9:]
            hgdp_cols = [9 + i for i, s in enumerate(samples) if s.startswith("HGDP")]
            hgdp_names = [samples[i - 9] for i in hgdp_cols]
            fout.write(
                "##INFO=<ID=AF_HGDP,Number=1,Type=Float,Description=\""
                "Alt allele frequency in HGDP samples in this VCF\">\n"
            )
            fout.write("##FORMAT=<ID=GT,Number=1,Type=String,"
                       "Description=\"Genotype\">\n")
            fout.write("\t".join(cols[:9] + hgdp_names) + "\n")
            continue

        # variant record
        n_in += 1
        f = line.rstrip("\n").split("\t")
        ref, alt = f[3], f[4]
        if len(ref) != 1 or len(alt) != 1 or "," in alt or alt == "*":
            continue
        n_pass += 1

        # HGDP-only genotypes (first ':'-delimited field of each col)
        ac = an = 0
        out_gts = []
        for i in hgdp_cols:
            g = f[i].split(":", 1)[0]
            if g in ("./.", ".|.", "."):
                out_gts.append("./.")
                continue
            a1, _, a2 = g.replace("|", "/").partition("/")
            if a1 in ("0", "1") and a2 in ("0", "1"):
                ac += int(a1) + int(a2)
                an += 2
                out_gts.append(f"{a1}/{a2}")
            else:
                out_gts.append("./.")
        if an == 0:
            continue
        af = ac / an
        if min(af, 1 - af) < MAF_MIN:
            continue
        n_kept += 1

        fout.write("\t".join([
            f[0], f[1], f[2], ref, alt, ".", "PASS",
            f"AF_HGDP={af:.5f}", "GT", *out_gts
        ]) + "\n")

        now = time.time()
        if now - last_print > PROGRESS_EVERY_S:
            print(f"  chr{chrom}: bytes_decompressed={bytes_seen/1e9:.2f} GB "
                  f"n_in={n_in:,} n_kept={n_kept:,}",
                  file=sys.stderr, flush=True)
            last_print = now

    return n_in, n_pass, n_kept


def stream_one(chrom: str) -> dict:
    dst = DDIR / f"chr{chrom}.hgdp_subset.vcf.gz"
    tmp = dst.with_suffix(dst.suffix + ".tmp")
    if dst.exists():
        return {"chrom": chrom, "action": "SKIP",
                "size_mb": round(dst.stat().st_size / 1e6, 1)}

    url = BASE + FILE.format(c=chrom)
    last_err = None
    for attempt in range(1, N_RETRIES + 1):
        if tmp.exists():
            tmp.unlink()
        t0 = time.time()
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "roh-calibration/0.1",
                    "Accept-Encoding": "identity",  # no transport gzip; we decode
                },
            )
            with urllib.request.urlopen(req, timeout=300) as resp, \
                 gzip.GzipFile(fileobj=resp, mode="rb") as gz, \
                 gzip.open(tmp, "wt", compresslevel=6) as fout:
                n_in, n_pass, n_kept = _filter_stream(gz, fout, chrom)
            tmp.rename(dst)
            return {
                "chrom": chrom, "action": "DONE",
                "size_mb": round(dst.stat().st_size / 1e6, 1),
                "n_in": n_in, "n_kept": n_kept,
                "elapsed_s": int(time.time() - t0),
                "attempt": attempt,
            }
        except (urllib.error.URLError, socket.timeout, ConnectionError,
                EOFError, OSError) as e:
            last_err = e
            wait = RETRY_BACKOFF_S * (2 ** (attempt - 1))
            print(f"chr{chrom}: attempt {attempt}/{N_RETRIES} failed "
                  f"({type(e).__name__}: {e}); retrying in {wait}s",
                  file=sys.stderr, flush=True)
            time.sleep(wait)

    return {"chrom": chrom, "action": "FAILED",
            "error": f"{type(last_err).__name__}: {last_err}"}


def main():
    chroms = sys.argv[1:] if len(sys.argv) > 1 else SMALLEST_FIRST
    print(f"streaming HGDP subset; chroms={chroms}", file=sys.stderr, flush=True)
    summary = []
    for c in chroms:
        print(f"\n=== chr{c} ===", file=sys.stderr, flush=True)
        r = stream_one(c)
        print(r, file=sys.stderr, flush=True)
        summary.append(r)

    print("\n=== summary ===", file=sys.stderr)
    for r in summary:
        print(r, file=sys.stderr)


if __name__ == "__main__":
    main()
