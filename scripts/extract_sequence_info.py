import csv
import sys
from collections import Counter

from Bio import SeqIO

reference_fasta_path = sys.argv[1]
output_path = sys.argv[2]


with open(reference_fasta_path) as reference_fasta, open(output_path, "w") as f:
    reference = SeqIO.parse(reference_fasta, "fasta")

    nucleotide_info = csv.writer(f)
    nucleotide_info.writerow(
        [
            "sequence_id",
            "a_percent",
            "c_percent",
            "t_percent",
            "gc_percent",
            "length",
        ]
    )
    
    for record in reference:
        counter = Counter(record.seq.lower())

        nucleotide_info.writerow(
            [   
                record.id,
                counter["a"] / counter.total(),
                counter["c"] / counter.total(),
                counter["t"] / counter.total(),
                (counter["g"] + counter["c"]) / counter.total(),
                len(record.seq),
            ]
        )
