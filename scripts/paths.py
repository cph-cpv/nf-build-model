from pathlib import Path
import sys


for stem in ("illumina", "rott"):
    for extension in ("fastq", "fq"):
        search_path = Path("/mnt/raw") / stem

        for p in search_path.rglob(f"*.{extension}*"):
            print(p)

