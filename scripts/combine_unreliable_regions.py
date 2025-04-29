import argparse
import csv

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine unreliable regions from two files.")

    parser.add_argument("nucleotide_regions", help="Path to the nucleotide regions file.")
    parser.add_argument("mapped_regions", help="Path to the mapped regions file.")
    parser.add_argument("output", help="Path to the output file.")

    args = parser.parse_args()

    with open(args.output, "w") as output_f:
        output = csv.writer(output_f)
        output.writerow(["sequence_id", "start", "end", "category"])

        with open(args.nucleotide_regions) as f:
            reader = csv.reader(f)
            next(reader) 
            output.writerows(([sequence_id,start,end,category] for start,end,sequence_id,category in reader))

        with open(args.mapped_regions) as f:
            reader = csv.reader(f)
            next(reader)

            output.writerows(
                [sequence_id, start, end, "Unmappable regions"] for sequence_id,start,end in reader
            )
