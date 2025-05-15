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

bam_files_path <- args[1]
unreliable_regions_path <- args[2]
nucleotide_info_path <- args[3]
sample_labels_path <- args[4]

nucleotide_info <- read.csv(nucleotide_info_path, header=TRUE)
unreliable_regions <- read.csv(unreliable_regions_path, header=TRUE, check.names=FALSE)

all_files_path <- list.files(bam_files_path)
rle_file_paths <- all_files_path[grep("^rle_.*.rds$", all_files_path)]

training_data <- lapply(rle_file_paths, function(rle_file_path){
        rle_data <- readRDS(file.path(bam_files_path, rle_file_path))
        convert_rle_to_df(rle_data, unreliable_regions, additional_nucleotide_info=nucleotide_info)
    })

training_data <- Filter(function(df) nrow(df) > 0, training_data)

train_x <- bind_rows(training_data)

sample_labels <- fromJSON(sample_labels_path)

train_y <- logical(nrow(train_x))

for(i in 1:nrow(train_x)){
    row <- train_x[i, ]
    train_y[i] <- row[["seg_id"]] %in% sample_labels[[row[["sample_id"]]]]
}

model <- train_iimi(train_x, train_y)

saveRDS(model, file = "trained_xgb.rds")
