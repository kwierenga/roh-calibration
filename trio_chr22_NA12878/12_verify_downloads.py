"""Stream-test gzip integrity of every downloaded chromosome VCF.
Reports the byte position where decompression fails (if any).
"""

import gzip
import sys
from pathlib import Path

BASE = Path(r"c:\Users\klaas\ROH calibration project\trio_chr22_NA12878\all_autosomes")
CHR22 = Path(r"c:\Users\klaas\ROH calibration project\trio_chr22_NA12878\chr22_phased.vcf.gz")

files = list(BASE.glob("chr*_phased.vcf.gz")) + [CHR22]
files.sort()

problems = []
for f in files:
    sz = f.stat().st_size
    decompressed = 0
    lines = 0
    try:
        with gzip.open(f, "rb") as fh:
            while True:
                chunk = fh.read(1024 * 1024)
                if not chunk:
                    break
                decompressed += len(chunk)
                lines += chunk.count(b"\n")
        print(f"  OK    {f.name:35s}  {sz/1e9:.2f} GB compressed -> {decompressed/1e9:.2f} GB raw, {lines:,} lines")
    except Exception as e:
        print(f"  FAIL  {f.name:35s}  failed at {decompressed:,} bytes in: {e}")
        problems.append(f.name)

if problems:
    print(f"\n  {len(problems)} corrupted file(s): {problems}")
else:
    print(f"\n  all {len(files)} files OK")
