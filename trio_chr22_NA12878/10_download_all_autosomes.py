"""
Download phased VCFs for chr1..chr21 from 1000G high-coverage 20220422 release,
in parallel (4 workers by default). chr22 already on disk and skipped.

Resumable: skips files already on disk at expected size.
No ENTER prompts; fully scripted; safe to run in background.
"""

import concurrent.futures as cf
import sys
import time
import urllib.request
from pathlib import Path

BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"

HERE = Path(__file__).parent
DATA_DIR = HERE / "all_autosomes"
DATA_DIR.mkdir(exist_ok=True)

CHROMS = [f"chr{n}" for n in range(1, 22)]  # chr22 already exists in HERE
PARALLEL = 4


def url_for(chrom):
    return f"{BASE}1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        return int(r.headers.get("Content-Length", "0"))


def download_one(chrom):
    url = url_for(chrom)
    dest = DATA_DIR / f"{chrom}_phased.vcf.gz"
    tbi_url = url + ".tbi"
    tbi_dest = DATA_DIR / f"{chrom}_phased.vcf.gz.tbi"

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

    # index
    try:
        if not tbi_dest.exists():
            urllib.request.urlretrieve(tbi_url, tbi_dest)
    except Exception:
        pass

    return chrom, action, f"{dest.stat().st_size:,} bytes"


def main():
    t0 = time.time()
    print(f"  downloading {len(CHROMS)} autosomes to {DATA_DIR}/")
    print(f"  parallel workers: {PARALLEL}")
    print()

    # sort by expected size descending so largest start first
    print("  inspecting sizes ...")
    sizes = {}
    for c in CHROMS:
        try:
            sizes[c] = remote_size(url_for(c))
        except Exception:
            sizes[c] = 0
    total_mb = sum(sizes.values()) / 1e6
    print(f"  total to download: {total_mb:,.0f} MB across {len(CHROMS)} files")
    chroms_sorted = sorted(CHROMS, key=lambda c: -sizes.get(c, 0))

    completed = 0
    with cf.ThreadPoolExecutor(max_workers=PARALLEL) as executor:
        futures = {executor.submit(download_one, c): c for c in chroms_sorted}
        for fut in cf.as_completed(futures):
            chrom, status, info = fut.result()
            completed += 1
            elapsed = time.time() - t0
            print(f"  [{completed:2d}/{len(CHROMS)}] {chrom:8s} {status:10s} {info}  ({elapsed:.0f}s elapsed)")
            sys.stdout.flush()

    elapsed = time.time() - t0
    print()
    print(f"  done in {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
