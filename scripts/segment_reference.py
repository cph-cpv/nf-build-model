import json
import copy
import sys 
"""
Takes a reference.json and write a fasta file of 75 bp seqs
"""

def write_sequence_fragments(sequence):
    id, sequence_nucleotides = sequence["_id"], sequence["sequence"]
    with open("scripts/output/fragmented.fasta", "w") as output:
        for start in range(0, len(sequence_nucleotides) - 75):
            output.write(f"> {id} \n")
            output.write(f"{sequence_nucleotides[start:start + 75]} \n\n")

with open("test_folder/reference.json","r") as referenceFile:
    reference = json.load(referenceFile)

    for otu in reference["otus"]:
        for isolate in otu["isolates"]:
            for sequence in isolate["sequences"]:
                write_sequence_fragments(sequence)
                





#     json.dump(reduced_reference, output)