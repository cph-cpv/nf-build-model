import csv
from typing import Iterator
from pprint import pprint
from pathlib import Path
import sys
import json

reference_path = sys.argv[1]


def parse_clstr(path: Path) -> Iterator[dict]:
    """Parse the clstr output file from cd-hit."""
    cluster = None

    with open(path) as f:
        for line in f:
            if line[0] == ">":
                if cluster:
                    yield cluster

                cluster = {
                    "id": line[1:].strip(),
                    "members": [],
                }

            else:
                sequence_id = line.split(">")[1].split("...")[0]
                cluster["members"].append(sequence_id)


def parse_fasta(path: Path):
    count = 0
    headers = set()

    with open(path) as f:
        for line in f:
            if line[0] == ">":
                headers.add(line[1:].strip())
                count += 1

    if len(headers) != count:
        raise ValueError("Duplicate headers found in FASTA file")

    return headers


def parse_clusters(path: Path):
    _reps_by_sequence_id = {}

    for output_path in path.iterdir():
        if not output_path.stem.startswith("output_"):
            continue

        representative_headers = parse_fasta(output_path / "clustered.fa")

        for cluster in parse_clstr(output_path / "clustered.fa.clstr"):
            reps = [s for s in cluster["members"] if s in representative_headers]

            if len(reps) != 1:
                raise ValueError(
                    f"Expected one representative for cluster {cluster['id']}"
                )

            for sequence_id in cluster["members"]:
                _reps_by_sequence_id[sequence_id] = reps[0]

    return _reps_by_sequence_id


reps_by_sequence_id = parse_clusters(Path.cwd())


all_edges = set()

clustered_isolate_count = 0
reference_isolate_count = 0

reference_segment_count = 0

clustered_sequence_count = len(reps_by_sequence_id)
reference_sequence_count = 0
representative_sequence_count = len(set(reps_by_sequence_id.values()))

reference_minimum_sequence_length = -1


with open(reference_path) as f:
    reference = json.load(f)


for otu in reference["otus"]:
    otu_id = otu["_id"]

    collapsed_isolates = []

    sequences_by_id = {
        sequence["_id"]: sequence
        for isolate in otu["isolates"]
        for sequence in isolate["sequences"]
    }

    reference_segment_count += len(otu["schema"])
    required_segments = {s["name"] for s in otu["schema"] if s["required"]}

    for isolate in otu["isolates"]:
        reference_isolate_count += 1

        reference_sequence_count += len(isolate["sequences"])

        reference_minimum_sequence_length = min(
            reference_minimum_sequence_length,
            min(len(sequence["sequence"]) for sequence in isolate["sequences"]),
        )

        edges = tuple(
            reps_by_sequence_id.get(sequence["_id"], sequence["_id"])
            for sequence in isolate["sequences"]
        )

        if edges in all_edges:
            continue

        all_edges.add(edges)

        sequences_by_segment = {
            sequences_by_id[sequence_id]["segment"]: sequences_by_id[sequence_id][
                "sequence"
            ]
            for sequence_id in edges
        }

        for sequence in isolate["sequences"]:
            segment_name = sequence["segment"]

            try:
                sequence["sequence"] = sequences_by_segment[segment_name]
            except KeyError as e:
                if segment_name in required_segments:
                    raise KeyError(
                        f"Missing required segment {segment_name} in sequences_by_segment for isolate {isolate['id']}"
                    ) from e

        collapsed_isolates.append(isolate)

    clustered_isolate_count += len(collapsed_isolates)

    otu["isolates"] = collapsed_isolates


with open("reps_by_sequence.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(["sequence_id", "representative_id"])

    for sequence_id, representative_id in reps_by_sequence_id.items():
        writer.writerow([sequence_id, representative_id])


written_rep_ids = set()


with open("reps.fa", "w") as f:
    for otu in reference["otus"]:
        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                sequence_id = sequence["_id"]

                representative_id = reps_by_sequence_id.get(sequence_id)

                if representative_id in written_rep_ids:
                    continue

                written_rep_ids.add(representative_id)

                f.write(f">{representative_id}\n{sequence['sequence']}\n")


with open("collapsed.json", "w") as f:
    json.dump(reference, f)


with open("summary.txt", "w") as f:
    f.write("Summary of clustering results:\n")
    f.write("==============================\n\n")

    f.write("Isolates:\n")
    f.write("---------\n\n")

    f.write(f"Input:     {reference_isolate_count}\n")
    f.write(f"Collapsed: {clustered_isolate_count}\n\n")

    f.write("Sequences:\n")
    f.write("----------\n\n")

    f.write(f"Input:     {reference_sequence_count}\n")
    f.write(f"Collapsed: {clustered_sequence_count}\n")
    f.write(f"Reps:      {representative_sequence_count}\n\n")

    f.write("Other:\n")
    f.write("------\n\n")

    f.write(f"Total OTU count:         {len(reference['otus'])}\n")
    f.write(f"Total segment count:     {reference_segment_count}\n")
    f.write(f"Minimum sequence length: {reference_minimum_sequence_length}\n\n")
