import csv
import json
import sys
from pathlib import Path
from pprint import pprint
from types import SimpleNamespace
from typing import Iterator

reference_path = sys.argv[1]


class CollapseCounts(SimpleNamespace):
    collapsed_isolates = 0
    reference_isolates = 0
    collapsed_sequences = 0
    reference_sequences = 0
    representative_sequences = 0


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
                if cluster is None:
                    raise ValueError("Unexpected line in clstr file")
                
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


def write_reps_by_sequence_id(reps_by_sequence_id_: dict[str, str]):
    with open("reps_by_sequence.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence_id", "representative_id"])

        for sequence_id, representative_id in reps_by_sequence_id_.items():
            writer.writerow([sequence_id, representative_id])


def write_reps(reference_: dict, reps_by_sequence_id_: dict[str, str]):
    written_sequence_ids = set()

    with open("reps.fa", "w") as f:
        for otu in reference_["otus"]:
            for isolate in otu["isolates"]:
                for sequence in isolate["sequences"]:
                    sequence_id = sequence["_id"]

                    if (
                        sequence_id in written_sequence_ids
                        or sequence_id in reps_by_sequence_id_
                    ):
                        continue

                    f.write(f">{sequence_id}\n{sequence['sequence']}\n")
                    written_sequence_ids.add(sequence_id)


def load_reference(
    path: Path, counts: CollapseCounts, reps_by_sequence_id_: dict[str, str]
) -> dict:
    with open(path) as f:
        reference = json.load(f)

    all_edges = set()
    reference_minimum_sequence_length = -1

    for otu in reference["otus"]:
        collapsed_isolates = []

        for isolate in otu["isolates"]:
            counts.reference_isolates += 1
            counts.reference_sequences += len(isolate["sequences"])

            reference_minimum_sequence_length = min(
                reference_minimum_sequence_length,
                min(len(sequence["sequence"]) for sequence in isolate["sequences"]),
            )

            edges = tuple(
                reps_by_sequence_id_.get(sequence["_id"], sequence["_id"])
                for sequence in isolate["sequences"]
            )

            if edges in all_edges:
                continue

            all_edges.add(edges)
            collapsed_isolates.append(isolate)

        counts.collapsed_isolates += len(collapsed_isolates)
        otu["isolates"] = collapsed_isolates

    return reference


if __name__ == "__main__":
    counts = CollapseCounts()

    reps_by_sequence_id = parse_clusters(Path.cwd())

    counts.collapsed_sequences = len(reps_by_sequence_id)
    counts.representative_sequences = len(set(reps_by_sequence_id.values()))

    reference = load_reference(Path(reference_path), counts, reps_by_sequence_id)
    write_reps_by_sequence_id(reps_by_sequence_id)
    write_reps(reference, reps_by_sequence_id)

    with open("collapsed.json", "w") as f:
        json.dump(reference, f)

    with open("summary.txt", "w") as f:
        f.write("Summary of clustering results:\n")
        f.write("==============================\n\n")

        f.write("Isolates:\n")
        f.write("---------\n\n")

        f.write(f"Input:     {counts.reference_isolates}\n")
        f.write(f"Collapsed: {counts.collapsed_isolates}\n\n")

        f.write("Sequences:\n")
        f.write("----------\n\n")

        f.write(f"Input:     {counts.reference_sequences}\n")
        f.write(f"Collapsed: {counts.collapsed_sequences}\n")
        f.write(f"Reps:      {counts.representative_sequences}\n\n")

        f.write("Other:\n")
        f.write("------\n\n")

        f.write(f"Total OTU count:         {len(reference['otus'])}\n")
        # f.write(f"Minimum sequence length: {reference_minimum_sequence_length}\n\n")
