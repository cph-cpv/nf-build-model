"""The the mapped regions in a BAM file."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pysam


def find_mapped_regions(bam_path: Path, output_path: Path, minimum_coverage: int = 1):
    """Find mapped regions in BAM files.

    Mapped regions are contiguous stretches of the reference genome that are mapped by at least one read.

    :param bam_path: the path to the BAM file
    """
    covered_positions: defaultdict[str, set[int]] = defaultdict(set)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam_file:
        for chromosome in bam_file.references:
            for read in bam_file.fetch(chromosome):
                if read.is_unmapped:
                    continue

                if read.query_name is None:
                    raise ValueError("Read has no query name")

                split = read.query_name.split(":")
                sequence_id = split[0]
                window_start = int(split[1])

                covered_positions[sequence_id].update(
                    range(
                        window_start + read.query_alignment_start - 1,
                        window_start + read.query_alignment_end,
                    )
                )

    with open(output_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence_id", "start", "end"])

        for sequence_id in covered_positions:
            positions = sorted(covered_positions[sequence_id])
            start = positions[0]
            end = positions[0]

            for pos in positions[1:]:
                if pos == end + 1:
                    end = pos
                else:
                    writer.writerow([sequence_id, start, end])
                    start = pos
                    end = pos


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find mapped regions in BAM files.")
    parser.add_argument("bam_path", type=Path, help="The path to the BAM file.")
    parser.add_argument(
        "output_path", type=Path, help="The path to the output CSV file."
    )
    parser.add_argument(
        "--minimum_coverage",
        type=int,
        default=1,
        help="Minimum coverage to consider a region mapped.",
    )

    args = parser.parse_args()
    find_mapped_regions(args.bam_path, args.output_path, args.minimum_coverage)
