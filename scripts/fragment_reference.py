import sys

from Bio import SeqIO

input_path = sys.argv[1]
output_path = sys.argv[2]


def _yield_fragments():
    for record in SeqIO.parse(input_path, "fasta"):
        print(record.id)

        if record.id == "None":
            raise ValueError("Found None record")

        for start in range(0, len(record.seq) - 75):
            yield SeqIO.SeqRecord(
                record.seq[start : start + 75],
                id=f"{record.id}:{start}",
                description="",
            )


with open(output_path, "w") as f:
    SeqIO.write(
        _yield_fragments(),
        f,
        "fasta",
    )
