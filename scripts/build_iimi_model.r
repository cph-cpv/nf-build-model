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

# fileNames <- args[1]
bam_files_path <- args[1]
unreliable_regions_path <- args[2]
nucleotide_info_path <- args[3]
sample_labels_path <- args[4]



nucleotide_info <- read.csv(nucleotide_info_path, header=TRUE)
unreliable_regions <- read.csv(unreliable_regions_path, header=TRUE)

sample_file_paths <- list.files(bam_files_path)
bam_file_paths <- sample_file_paths[grep(".bam$", sample_file_paths)]


training_data <- lapply(bam_file_paths, function(bam_file_path){
        rle_data <- convert_bam_to_rle(file.path(bam_files_path, bam_file_path))
        convert_rle_to_df(rle_data, unreliable_regions, additional_nucleotide_info=nucleotide_info)
    })

train_x <- bind_rows(training_data)


sample_labels <- fromJSON(sample_labels_path)

train_y <- logical(nrow(train_x))

for(i in 1:nrow(train_x)){
    row <- train_x[i, ]
    train_y[i] <- row[["seg_id"]] %in% sample_labels[[row[["sample_id"]]]]
}


model <- train_iimi(train_x, train_y)

saveRDS(model, file = "trained_xgb.rds")
