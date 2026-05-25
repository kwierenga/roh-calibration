"""
Download two external datasets needed for the next methodology iteration:

1. Palsson et al. 2024 Nature deCODE "Complete human recombination maps"
   (a successor to Halldorsson 2019 with full coverage and improved resolution).
   Sex-specific maps in 8.1 MB zip.

2. Illumina Platinum Genomes 2017-1.0 release: NA12878 hg38 phased VCF.
   Pedigree-phased against the full 17-member CEPH 1463 pedigree.
   Will let us compare 1000G-statistical-phased NA12878 vs Platinum-pedigree-
   phased NA12878 on chr22 to localize where statistical phasing disagrees
   with pedigree truth.

Both small. Resumable if already on disk at expected size.
"""

import sys
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).parent
EXTERNAL_DIR = HERE / "external"
EXTERNAL_DIR.mkdir(exist_ok=True)

FILES = {
    "palsson2024_deCODE_maps.zip": "https://zenodo.org/records/14025565/files/DecodeGenetics/PalssonEtAl_Nature_2024-v1.0.0.zip",
    "platinum_NA12878_hg38.vcf.gz": "https://platinum-genomes.s3.amazonaws.com/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz",
}


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", "0"))


def download(name, url):
    dest = EXTERNAL_DIR / name
    try:
        expected = remote_size(url)
    except Exception as e:
        print(f"  [warn] HEAD failed for {name}: {e}")
        expected = -1

    if dest.exists() and expected > 0 and dest.stat().st_size == expected:
        print(f"  [skip] {name} already at {expected:,} bytes")
        return dest

    print(f"  [get ] {name} ({expected/1e6:.1f} MB) <- {url}")

    def hook(blocks, blocksize, total):
        if total <= 0:
            return
        downloaded = min(blocks * blocksize, total)
        pct = 100 * downloaded / total
        bar = "#" * int(pct / 2.5)
        sys.stdout.write(f"\r    [{bar:<40}] {pct:5.1f}%")
        sys.stdout.flush()

    urllib.request.urlretrieve(url, dest, reporthook=hook)
    print()
    return dest


def main():
    paths = {}
    for name, url in FILES.items():
        paths[name] = download(name, url)

    # unzip the deCODE archive
    decode_zip = paths["palsson2024_deCODE_maps.zip"]
    decode_dir = EXTERNAL_DIR / "palsson2024_deCODE_maps"
    if decode_zip.exists() and not decode_dir.exists():
        print(f"  unzipping {decode_zip.name} -> {decode_dir}")
        with zipfile.ZipFile(decode_zip) as zf:
            zf.extractall(decode_dir)
        print(f"  contents:")
        for p in sorted(decode_dir.rglob("*")):
            if p.is_file():
                print(f"    {p.relative_to(decode_dir)}  ({p.stat().st_size:,} bytes)")

    print("done.")


if __name__ == "__main__":
    main()
