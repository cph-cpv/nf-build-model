"""Scan a given path for files with specific extensions."""

from pathlib import Path

for stem in ("illumina", "rott"):
    for extension in ("fastq", "fq"):
        search_path = Path("/mnt/raw") / stem

        for p in search_path.rglob(f"*.{extension}*"):
            print(p)

