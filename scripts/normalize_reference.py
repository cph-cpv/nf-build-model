import csv
import json
from collections import Counter, defaultdict

sample_data = defaultdict(set)

with open("data/cleaned.csv") as f:
    reader = csv.reader(f)    

    for row in reader:
        sample_id, extract, viruses = row

        for acronym in viruses.split(","):
            sample_data[sample_id].add(acronym)

sorted_sample_ids = sorted(sample_data.keys())

with open("test_folder/reference.json") as f:
    reference_json = json.load(f)


seen_acronyms = set()


with open("scripts/output/nucleotide_info.csv", "w") as f_ni, open("scripts/output/labels.csv", "w") as f_l:
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

    labels = csv.writer(f_l)
    labels.writerow(["sequence_ids", *sorted_sample_ids])

    abbr_to_otu_ids = defaultdict(list)

    for otu in reference_json["otus"]:
        otu_id = otu["_id"]

        acronym = otu["abbreviation"].lower()
        name = otu["name"].lower()

        for char in ("_", "-", " "):
            acronym.replace(char, "")

        if name == "peach chlorotic mottle virus":
            acronym = "pecmv"

        if name == "cucumber fruit mottle mosaic virus":
            acronym = "cufmmv"

        if acronym:
            abbr_to_otu_ids[acronym].append(otu_id)

        if acronym and acronym in seen_acronyms:
            raise ValueError(f"Already saw acronym: {acronym}")

        seen_acronyms.add(acronym)
        
        for isolate in otu["isolates"]:
            isolate_id = isolate["id"]

            for sequence in isolate["sequences"]:
                sequence_id = sequence["_id"]

                counter = Counter(sequence["sequence"].lower())

                a_percent = counter["a"] / counter.total()

                nucleotide_info.writerow(
                    [
                        otu_id,
                        acronym,
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

                row = [
                    acronym in sample_data[sample_id] for sample_id in sorted_sample_ids
                ]

                labels.writerow([sequence_id, *row])

                if all(value is False for value in row):
                    print("wrote falsey row")
