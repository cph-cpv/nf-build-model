library(Biostrings)
library(iimi)
library(mltools)
library(data.table)
library(dplyr)
library(randomForest)
library(Rsamtools)
library(GenomicAlignments)
library(jsonlite)

args <- commandArgs(trailingOnly=TRUE)

fileName <- args[1]
reference_json_path <- args[2]
virus_segments_path <-  args[3]
host_mapping_path <- args[4]
nucleotide_info_path <- args[5]

nucleotide_info <- read.csv(nucleotide_info_path, sep="\t", header=TRUE)

# reference = fromJSON(reference_json_path)

virus_segments = readDNAStringSet(virus_segments_path)

# head(virus_segments)

#Create unreliable regions
unreliable_regions = create_mappability_profile(host_mapping_path, "Unmappable region (Host)", 75, virus_segments )
#  head(unreliable_regions)

rle_data = convert_bam_to_rle(paste(fileName, ".bam", sep=""))



for(row in rle_data){
    row
    seq_id = row$
}
 

save(object=rle_data, file=paste(fileName, "_rle_data.rds", sep=""))

# bam_file_paths <- file_paths[grep(".bam$", file_paths)]
# bam_file_paths <- lapply(bam_file_paths, function(file_name){
#     file.path("./", bam_path, file_name)
# })

# bam_file_paths

# rle_data <- lapply( bam_file_paths, convert_bam_to_rle)





