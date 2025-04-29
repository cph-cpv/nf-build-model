import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

EXCLUDED_SAMPLE_NAMES = [
    "QUADS39-rep2-smRNA_GTAGCC_R2_001.fq.gz",
"QUADS40-rep1-smRNA_ATTGGC_R1_001.fq.gz",
"QUADS40-rep1-smRNA_ATTGGC_R2_001.fq.gz",
]


def check_labels_without_fastq(sample_names: set[str], seen: dict[str, set[Path]]) -> None:
    if labels_without_fastq := [name for name in sample_names if name not in seen]:
        joined = ", ".join([f"'{name}'" for name in sorted(labels_without_fastq)])
        print(
            "Found sample labels without corresponding FASTQ files. Make sure all labels have "
            f"corresponding FASTQ files: {joined}",
            file=sys.stderr, 
        )

        sys.exit(1)


def associate(labels: list[dict[str, str]], sample_paths: list[Path]) -> None:
    sample_names = {label["sample_name"] for label in labels if label["sample_name"] not in EXCLUDED_SAMPLE_NAMES}

    seen: dict[str, set[Path]] = defaultdict(set)


    for sample_name in sample_names:
        for path in sample_paths:
            if sample_name in path.name:
                seen[sample_name].add(path)

    # check_labels_without_fastq(sample_names, seen)

    if fastq_without_labels := (
        set(sample_paths) - {path for paths in seen.values() for path in paths}
    ):
        joined = ", ".join([f"'{path.name}'" for path in sorted(fastq_without_labels)])

        print(
            "Found FASTQ files without sample labels. Make sure all files have corresponding "
            f"labels: {joined}"
        )

        sys.exit(1)
    
    if non_unique_sample_names:=[
        (sample_name, paths) for sample_name, paths in seen.items() if len(paths) > 1
    ]:
        print(
            "Found non-specific sample names. Update your sample labels to be more specific to "
            "your sample FASTQ files.\n"
        )

        for sample_name, paths in non_unique_sample_names:
            if len(seen[sample_name]) > 1:
                print(f"{sample_name}:")

                for path in paths:
                    print(f"  • {path}")

        sys.exit(1)

    rows = []

    for label in labels:
        sample_name = label["sample_name"]

        for path in sample_paths:
            if sample_name in path.name:
                rows.append((label["virus_id"], sample_name, path))

    with open(args.output_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["virus_id", "sample_name", "path"])
        writer.writerows(rows)


def get_sample_paths(samples_dir_path: Path) -> list[Path]:
    sample_paths = []

    for pattern in ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]:
        for path in samples_dir_path.rglob(pattern):
            if path.is_dir():
                raise ValueError(f"Path is a directory: {path}")

            sample_paths.append(path)

    return sample_paths


def parse_sample_labels(path: Path) -> list[dict[str, str]]:
    with open(path) as f:
        return [
            {"virus_id": row["virus_id"], "sample_name": row["sample_name"]}
            for row in csv.DictReader(f) if row["sample_name"] not in EXCLUDED_SAMPLE_NAMES
        ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a CSV file mapping virus call labels to files in a sample directory."
    )
    parser.add_argument(
        "labels_path",
        type=Path,
        help="The path to the CSV file containing virus call labels.",
    )
    parser.add_argument(
        "samples_dir_path",
        type=Path,
        help="The path to directory containing FASTQ files.",
    )
    parser.add_argument(
        "output_path", type=Path, help="The path to the output CSV file."
    )

    args = parser.parse_args()

    labels = parse_sample_labels(args.labels_path)
    sample_paths = get_sample_paths(args.samples_dir_path)
    
    associate(labels, sample_paths)
