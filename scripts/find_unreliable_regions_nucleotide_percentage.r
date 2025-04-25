library(Biostrings)
library(iimi)

args <- commandArgs(trailingOnly=TRUE)

virus_segments_path <-  args[1]
output_path <- args[2]


virus_segments = readDNAStringSet(virus_segments_path)

unreliable_regions = create_high_nucleotide_content(0.6, 0.45, 75, virus_segments)

write.csv(unreliable_regions, output_path, row.names=FALSE)