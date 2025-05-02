import json
import sys
import csv

sample_labels_path = sys.argv[1]
reference_json_path = sys.argv[2]
output_path = sys.argv[3]

with open (reference_json_path, "r") as f:
    reference = json.load(f)

    virus_segs = {} #{virus_id: [sequence_id_1, ...],}

    for otu in reference["otus"]:
        segments = []
        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                segments.append(sequence["_id"])
        virus_segs[otu["_id"]] = segments


with open(sample_labels_path, "r") as sample_labels_f, open(output_path, "w") as output_f:
    sample_labels = csv.reader(sample_labels_f)
    next(sample_labels)
    
    sample_viruses = {}
   
    for row in sample_labels:
        sample_viruses.setdefault(row[1], []).extend(virus_segs[row[0]])

    json.dump(sample_viruses, output_f)