import csv
import json


extractions = set()
cleaned = []
all_viruses = set()

with open("input/samples.csv") as f:
    reader = csv.reader(f)

    for row in reader:
        name, _, extraction, viruses = row

        if name.endswith("A/B/C"):
            for prefix in ("A", "B", "C"):
                cleaned.append((name.replace("A/B/C", prefix), extraction, viruses))

        else:
            cleaned.append((name, extraction, viruses))

        all_viruses.update(viruses.rstrip().split(","))
        extractions.add(extraction)

for virus in sorted(all_viruses):
    print(virus)

with open("cleaned.csv", "w") as f:
    writer = csv.writer(f, quoting=csv.QUOTE_ALL)

    for sample in cleaned:
        writer.writerow(sample)
