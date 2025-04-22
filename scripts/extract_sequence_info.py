import csv
import json
from collections import Counter, defaultdict
import sys

reference_path = sys.argv[1]

with open(reference_path) as f:
    reference_json = json.load(f)


with open("scripts/output/nucleotide_info.csv", "w") as f:
    nucleotide_info = csv.writer(f)
    nucleotide_info.writerow(
        [
            "sequence_id",
            "virus_acronym",
            "virus_id",
            "virus_name",
            "isolate_id",
            "a_percent",
            "c_percent",
            "t_percent",
            "gc_percent",
            "length",
        ]
    )
    
    for otu in reference_json["otus"]:

        for isolate in otu["isolates"]:
            isolate_id = isolate["id"]

            for sequence in isolate["sequences"]:
                counter = Counter(sequence["sequence"].lower())

                nucleotide_info.writerow(
                    [   
                        sequence["_id"],
                        otu["abbreviation"]
                        otu["_id"],
                        otu["name"].lower(),
                        isolate_id,
                        counter["a"] / counter.total(),
                        counter["c"] / counter.total(),
                        counter["t"] / counter.total(),
                        (counter["g"] + counter["c"]) / counter.total(),
                        len(sequence["sequence"]),
                    ]
                )
