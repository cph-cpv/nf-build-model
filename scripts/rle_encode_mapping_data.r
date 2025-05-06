library(iimi)
library(Rsamtools)
library(GenomicAlignments)

args <- commandArgs(trailingOnly=TRUE)

bam_path <- args[1]
output_path <- args[2]
sample_name <- sub(pattern = ".sorted.bam$", replacement = "", basename(bam_path))

rle_data <- convert_bam_to_rle(bam_path)

saveRDS(rle_data, file = file.path(output_path, paste0("rle_", sample_name, ".rds")))
