"""
Download chr22 phased VCF and pedigree file from 1000G high-coverage release
(Byrska-Bishop et al. 2022, Cell).

Resumable: skips files already on disk at expected size.
"""

import sys
import urllib.request
from pathlib import Path

BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
FILES = {
    "chr22_phased.vcf.gz": BASE + "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
    "chr22_phased.vcf.gz.tbi": BASE + "20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi",
    "pedigree.txt": BASE + "1kGP.3202_samples.pedigree_info.txt",
}

HERE = Path(__file__).parent


def remote_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", "0"))


def download(url: str, dest: Path) -> None:
    expected = remote_size(url)
    if dest.exists() and dest.stat().st_size == expected:
        print(f"  [skip] {dest.name} already at {expected:,} bytes")
        return

    print(f"  [get ] {dest.name} ({expected/1e6:.1f} MB) <- {url}")

    def hook(blocks, blocksize, total):
        if total <= 0:
            return
        downloaded = min(blocks * blocksize, total)
        pct = 100 * downloaded / total
        bar = "#" * int(pct / 2.5)
        sys.stdout.write(f"\r    [{bar:<40}] {pct:5.1f}%  {downloaded/1e6:.1f}/{total/1e6:.1f} MB")
        sys.stdout.flush()

    urllib.request.urlretrieve(url, dest, reporthook=hook)
    print()


def main():
    for name, url in FILES.items():
        download(url, HERE / name)
    print("done.")


if __name__ == "__main__":
    main()
