import csv
from pathlib import Path
from collections import defaultdict
import sys
import re
import os


names = set()


with open("sample_names.txt") as f:
    for name in f:
        name = name.rstrip()

        if name:
            if name in names:
                raise ValueError(f"Non-unique name: {name}")

            names.add(name)


with open("fastq_paths.txt") as f:
    fastq_paths = [path.strip() for path in f if path]


paths_by_sample = defaultdict(list)


for path in fastq_paths:
    lowered = path.lower()

    for name in names:
        if name.lower() in lowered:
            paths_by_sample[name].append(path)


def resolvePaths(paths):
    file = {}
    for path in paths:
        if re.search("R2(_001)?\.fastq\.gz$", path):
            pass
        
        size = os.path.getsize(path)

        if file.get("size", 0) < size:
            file["path"] = path
            file["size"] = size
    
    if not file.get("path"):
        raise ValueError("No valid files found")

    return [file["path"]]


#Clean the paths to remove duplicates
for sample in paths_by_sample:
    if len(paths_by_sample[sample]) > 1:
        paths_by_sample[sample] = resolvePaths(paths_by_sample[sample])


with open("sample_name_paths.csv", "w") as f:
    writer = csv.writer(f)


    for sample_name, paths in paths_by_sample.items():
        for path in paths:
            writer.writerow([sample_name, path])


