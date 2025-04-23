"""The the mapped regions in a BAM file."""

import argparse
import pysam
import csv
from pathlib import Path


def find_mapped_regions(bam_path: Path, output_path: Path, minimum_coverage: int = 1):
    """Find mapped regions in BAM files.

    Mapped regions are contiguous stretches of the reference genome that are mapped by at least one read.

    :param path: the path to the BAM file
    """
    mapped_regions = []

    with pysam.AlignmentFile(bam_path) as bam_file:
        for sequence_id in bam_file.references:
            start = None

            for column in bam_file.pileup(sequence_id, stepper="all"):
                if column.n >= minimum_coverage and start is None:
                    start = column.reference_pos

                elif column.n < minimum_coverage:
                    if start is not None:
                        mapped_regions.append(
                            (sequence_id, start, column.reference_pos - 1)
                        )
                        start = None

    with open(output_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence_id", "start", "end"])
        writer.writerows(mapped_regions)


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
