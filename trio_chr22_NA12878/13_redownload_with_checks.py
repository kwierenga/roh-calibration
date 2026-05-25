"""
Re-download specific chromosomes single-threaded with post-download
integrity check. Retries up to N times on corruption.

Parallel-4 downloads introduced gzip-stream corruption in chr3 and chr10
during the initial bulk run; this single-threaded approach is slower but
the per-file gzip stream stays intact.
"""

import gzip
import sys
import time
import urllib.request
from pathlib import Path

BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
HERE = Path(__file__).parent
DATA_DIR = HERE / "all_autosomes"

CHROMS_TO_FIX = ["chr3", "chr10"]
MAX_RETRIES = 3


def url_for(chrom):
    return f"{BASE}1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def verify_gzip(path):
    try:
        with gzip.open(path, "rb") as fh:
            while True:
                chunk = fh.read(4 * 1024 * 1024)
                if not chunk:
                    return True
    except Exception as e:
        print(f"    integrity FAIL: {e}")
        return False


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        return int(r.headers.get("Content-Length", "0"))


def download_with_verify(chrom):
    url = url_for(chrom)
    dest = DATA_DIR / f"{chrom}_phased.vcf.gz"
    expected = remote_size(url)
    print(f"  {chrom}: target size {expected/1e9:.2f} GB")
    sys.stdout.flush()

    for attempt in range(1, MAX_RETRIES + 1):
        t0 = time.time()
        print(f"    attempt {attempt}/{MAX_RETRIES}: downloading single-threaded...")
        sys.stdout.flush()
        try:
            urllib.request.urlretrieve(url, dest)
        except Exception as e:
            print(f"    download error: {e}")
            continue
        elapsed = time.time() - t0
        actual = dest.stat().st_size
        if actual != expected:
            print(f"    size mismatch: got {actual:,}, expected {expected:,}")
            continue
        print(f"    downloaded {actual:,} bytes in {elapsed:.0f}s; verifying integrity...")
        sys.stdout.flush()
        if verify_gzip(dest):
            print(f"    {chrom}: OK")
            return True
        print(f"    {chrom}: corruption detected; retrying")
    print(f"    {chrom}: FAILED after {MAX_RETRIES} attempts")
    return False


def main():
    t0 = time.time()
    results = {}
    for chrom in CHROMS_TO_FIX:
        results[chrom] = download_with_verify(chrom)
    print()
    print(f"  done in {time.time() - t0:.0f}s")
    for chrom, ok in results.items():
        print(f"    {chrom}: {'OK' if ok else 'FAILED'}")
    if all(results.values()):
        print("  ALL_OK")
    else:
        print("  PARTIAL")


if __name__ == "__main__":
    main()
