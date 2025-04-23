"""
Takes a reference.json file and writes a new reference.json
with only the default isolates

"""

import json
import copy


with open("reference.json", "r") as referenceFile:
    reference = json.load(referenceFile)
    reduced_reference = copy.deepcopy(reference)


for otu in reduced_reference["otus"]:
    for isolate in otu["isolates"]:
        if isolate["default"]:
            otu["isolates"] = [isolate]
            break


with open("reduced_reference.json", "w") as output:
    json.dump(reduced_reference, output)
