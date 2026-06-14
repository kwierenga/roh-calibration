"""
Download gnomAD v3.1.2 HGDP+1KG joint callset (per-autosome VCFs) for HGDP
populations relevant to consanguinity validation:

  MENA:      Bedouin (46), Druze (42), Palestinian (46), Mozabite (27)
  Pakistani: Pathan, Sindhi, Burusho, Hazara, Brahui, Makrani, Balochi, Kalash
             (~22-25 each, ~200 total)
  Plus:      retains all HGDP populations for SGDP-style geographic overlay

Source: gnomAD public bucket on Google Cloud Storage (HTTPS, no auth).
URL pattern:
  https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/
  genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chrN.vcf.bgz

Resumable: skips files already on disk at expected remote size.
"""

import concurrent.futures as cf
import sys
import time
import urllib.request
from pathlib import Path

BASE = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
        "release/3.1.2/vcf/genomes/")
FILE = "gnomad.genomes.v3.1.2.hgdp_tgp.chr{c}.vcf.bgz"
TBI = FILE + ".tbi"

HERE = Path(__file__).parent
DATA_DIR = HERE / "hgdp_tgp"
DATA_DIR.mkdir(exist_ok=True)

CHROMS = sys.argv[1:] if len(sys.argv) > 1 else [str(n) for n in range(1, 23)]
PARALLEL = 3   # gnomAD bucket throughput; 3 workers is conservative


def url_for(chrom, ext="vcf"):
    f = FILE if ext == "vcf" else TBI
    return BASE + f.format(c=chrom)


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        return int(r.headers.get("Content-Length", "0"))


def download_one(chrom):
    url = url_for(chrom)
    dest = DATA_DIR / f"chr{chrom}.hgdp_tgp.vcf.bgz"
    tbi_url = url_for(chrom, "tbi")
    tbi_dest = DATA_DIR / f"chr{chrom}.hgdp_tgp.vcf.bgz.tbi"

    try:
        expected = remote_size(url)
    except Exception as e:
        return chrom, "ERROR_HEAD", str(e)

    if dest.exists() and dest.stat().st_size == expected:
        action = "SKIP"
    else:
        try:
            urllib.request.urlretrieve(url, dest)
            action = "OK"
        except Exception as e:
            return chrom, "ERROR_GET", str(e)

    try:
        if not tbi_dest.exists():
            urllib.request.urlretrieve(tbi_url, tbi_dest)
    except Exception:
        pass

    return chrom, action, f"{dest.stat().st_size:,} bytes"


def main():
    t0 = time.time()
    print(f"  HGDP+1KG joint callset -> {DATA_DIR}/")
    print(f"  chroms: {','.join('chr'+c for c in CHROMS)}")
    print(f"  parallel workers: {PARALLEL}")
    print()

    print("  inspecting remote sizes ...")
    sizes = {}
    for c in CHROMS:
        try:
            sizes[c] = remote_size(url_for(c))
            print(f"    chr{c}: {sizes[c]/1e6:,.0f} MB")
        except Exception as e:
            sizes[c] = 0
            print(f"    chr{c}: ERROR ({e})")
    total_mb = sum(sizes.values()) / 1e6
    print(f"  total: {total_mb:,.0f} MB across {len(CHROMS)} files")
    print()

    order = sorted(CHROMS, key=lambda c: -sizes.get(c, 0))
    with cf.ThreadPoolExecutor(max_workers=PARALLEL) as ex:
        for chrom, action, msg in ex.map(download_one, order):
            print(f"  [chr{chrom}] {action}: {msg}  ({time.time()-t0:.0f}s)")
            sys.stdout.flush()
    print(f"\n  total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
