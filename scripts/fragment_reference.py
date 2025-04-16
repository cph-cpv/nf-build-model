import sys

from Bio import SeqIO, SeqRecord

input_path = sys.argv[1]
output_path = sys.argv[2]



with open(output_path, "w") as f:
    def _yield_fragments():
        for record in SeqIO.parse(input_path, "fasta"):
            for start in range(0, len(record.seq) - 75):
                yield SeqRecord(record.seq[start:start + 75], id=f"{record.id}:{start}")

    SeqIO.write(
        _yield_fragments(),
        f,
        "fasta",
    )
