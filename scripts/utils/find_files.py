import csv
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def find_fastq_paths(path: Path, stems: list[str]) -> list[Path]:
    cache_path = Path.cwd() / ".cache"
    cache_path.mkdir(parents=True, exist_ok=True)

    try:
        with open(cache_path / "fastq_paths.json", "r") as f:
            cache = json.load(f)
            return [Path(p) for p in cache["paths"]]
    except FileNotFoundError:
        pass


    paths = []

    for stem in stems:
        for extension in ("fastq", "fq"):
            for p in (path / stem).rglob(f"*.{extension}*"):
                paths.append(p)
    
    with open(cache_path / "fastq_paths.json", "w") as f:
        json.dump({"paths": [str(path) for path in paths]}, f)

    return paths


def parse_sample_names(path: Path) -> set[str]:
    with open(path) as f:
        names = {row["sample_name"] for row in csv.DictReader(f)}

    if len(names) != len(set(names)):
        print("Warning: Non-unique sample names found.\n")
        sys.exit(1)

    return names


if __name__ == "__main__":
    names = parse_sample_names(Path.cwd() / "input/sample_labels.csv")
    paths = find_fastq_paths(Path("/mnt/raw"), ["illumina", "rott"])

    for path in (Path.cwd() / "input" / "samples").iterdir():
        if path.is_symlink():
            path.unlink()
    
    for name in names:        
        specific_paths = sorted([
            path
            for path in paths
            if name in path.name and "_R2" not in path.name
        ])

        try:

            path = specific_paths[0]
        except IndexError:
            print(f"Warning: No path found for sample name {name}.")
            continue

        (Path.cwd() / "input" / "samples" / path.name).symlink_to(path)
