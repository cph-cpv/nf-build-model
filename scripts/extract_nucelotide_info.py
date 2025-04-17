import csv
import json
from collections import Counter, defaultdict
import sys

reference_path = sys.argv[1]

with open(reference_path) as f:
    reference_json = json.load(f)


with open("scripts/output/nucleotide_info.csv", "w") as f_ni:
    nucleotide_info = csv.writer(f_ni)
    nucleotide_info.writerow(
        [
            "virus_name",
            "abbr",
            "name",
            "iso_id",
            "seg_id",
            "A_percent",
            "C_percent",
            "T_percent",
            "GC_percent",
            "seg_len",
        ]
    )

    for otu in reference_json["otus"]:
        otu_id = otu["_id"]
        name = otu["name"].lower()

        for isolate in otu["isolates"]:
            isolate_id = isolate["id"]

            for sequence in isolate["sequences"]:
                sequence_id = sequence["_id"]

                counter = Counter(sequence["sequence"].lower())

                a_percent = counter["a"] / counter.total()

                nucleotide_info.writerow(
                    [
                        otu_id,
                        name,
                        isolate_id,
                        sequence_id,
                        a_percent,
                        counter["c"] / counter.total(),
                        counter["t"] / counter.total(),
                        (counter["g"] + counter["c"]) / counter.total(),
                        len(sequence["sequence"]),
                    ]
                )
