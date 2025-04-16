import json


def correct_by_defintion(otu, corrections: list[tuple[str, str, str]]):
    """
    Correct the segment of sequences in the given OTU based on the provided corrections.
    """
    for name, in_definition, replacement in corrections:
        if otu["name"] == name:
            for isolate in otu["isolates"]:
                for sequence in isolate["sequences"]:
                    if in_definition in sequence["definition"]:
                        sequence["segment"] = replacement
                        break


with open("reference.json") as f:
    reference = json.load(f)

for otu in reference["otus"]:
    correct_by_defintion(
        otu,
        [
            (
                "Broad bean stain virus",
                "RNA2",
                "RNA 2",
            ),
            (
                "Actinidia chlorotic ringspot associated virus",
                "RNA3",
                "RNA 3",
            ),
            (
                "Citrus chlorotic spot virus",
                "RNA2",
                "RNA 2",
            ),
            (
                "Peach leaf pitting-associated virus",
                "RNA2",
                "RNA 2",
            ),
            (
                "Tomato golden vein virus",
                "DNA-B",
                "DNA B",
            ),
        ],
    )


with open("reference_repaired.json", "w") as f:
    json.dump(reference, f)
