import sys
import csv


with open(sys.argv[1]) as nucleotide_regions_f,  open(sys.argv[2], "r") as mapped_regions_f:
    
    mapped_regions = csv.reader(mapped_regions_f)
    nucleotide_regions = csv.reader(nucleotide_regions_f)
    next(mapped_regions), next(nucleotide_regions)

    with open(sys.argv[3], "w") as output_f:
        output = csv.writer(output_f)
        output.writerow(["sequence_id", "start", "end"])
        
        for region in mapped_regions:
            output.writerow(region)
        
        for region in nucleotide_regions:
            output.writerow(region[0:3])