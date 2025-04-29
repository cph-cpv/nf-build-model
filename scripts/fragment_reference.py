import argparse
from pathlib import Path
from typing import Iterator

from Bio import SeqIO


def create_fragments(input_path: Path) -> Iterator[SeqIO.SeqRecord]:
    for record in SeqIO.parse(input_path, "fasta"):
        if record.id == "None":
            raise ValueError("Found None record")

        for start in range(0, len(record.seq) - 75):
            yield SeqIO.SeqRecord(
                record.seq[start : start + 75],
                id=f"{record.id}:{start}",
                description="",
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fragment a reference genome into 75bp fragments"
    )
    parser.add_argument(
        "input_path",
        type=str,
        help="Path to the input reference genome in FASTA format",
    )
    parser.add_argument(
        "output_path",
        type=str,
        help="Path to the output file to save the fragments",
    )
    args = parser.parse_args()

    with open(args.output_path, "w") as f:
        SeqIO.write(
            create_fragments(args.input_path),
            f,
            "fasta",
        )
