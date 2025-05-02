import csv
import json
from collections import Counter, defaultdict
import sys
from Bio import SeqIO
import json

sequence_info_path = sys.argv[1]
reference_json_path = sys.argv[2]
output_path = sys.argv[3]

with open (reference_json_path, "r") as f:
    reference = json.load(f)

    seq_info = {} # seq_id: {sequence_id: {virus_name, isolate_id},}

    for otu in reference["otus"]:
        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                seq_info[sequence["_id"]] = {"virus_name": otu["name"], "isolate_id": isolate["id"]}


with open(sequence_info_path) as sequence_info_f, open(output_path, "w") as output_f:
    sequence_info = csv.reader(sequence_info_f)
    next(sequence_info)

    output = csv.writer(output_f)
    output.writerow(
        [
            "virus_name",
            "iso_id",
            "seg_id",
            "a_percent",
            "c_percent",
            "t_percent",
            "gc_percent",
            "length",
        ]
    )

    
    for row in sequence_info:
        seg_id = row.pop(0)

        output.writerow(
            [   seq_info[seg_id]["virus_name"],
                seq_info[seg_id]["isolate_id"],
                seg_id,
                *row
            ]
        )
