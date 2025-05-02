import json
import sys
import csv

sample_labels_path = sys.argv[1]
output_path = sys.argv[2]

with open(sample_labels_path, "r") as sample_labels_f, open(output_path, "w") as output_f:
    sample_labels = csv.reader(sample_labels_f)
    next(sample_labels)
    
    sample_viruses = {}
   
    for row in sample_labels:
        sample_viruses.setdefault(row[1], []).append(row[0])

    json.dump(sample_viruses, output_f)