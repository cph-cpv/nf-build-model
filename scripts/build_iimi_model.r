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

#import nucleotide info
#import unreliable regions
#do conversion for each bamfile

nucleotide_info <- read.csv(nucleotide_info_path, header=TRUE)
# nucleotide_info$v_name <- NULL
# nucleotide_info$isolate_id <- NULL

# head(nucleotide_info)

# nucleotide_info <- nucleotide_info[, c(7,8,1,2,3,4,5,6)]



unreliable_regions <- read.csv(unreliable_regions_path, header=TRUE)

sample_file_paths <- list.files(bam_files_path)

bam_file_paths <- sample_file_paths[grep(".bam$", sample_file_paths)]






training_data = lapply(bam_file_paths, function(bam_file_path){
        rle_data = convert_bam_to_rle(file.path(bam_files_path, bam_file_path))

        for (sample in names(rle_data)) {
            for (seg in names(rle_data[[sample]])) {
            isolate_id <- nucleotide_info[nucleotide_info$seg_id == seg, 2]
            v_name <- nucleotide_info[nucleotide_info$seg_id == seg, 1]
            print(paste("seg? :",seg, isolate_id))
            }
        }

        df_data = convert_rle_to_df(rle_data, unreliable_regions, additional_nucleotide_info=nucleotide_info)

        return(df_data)
    })

sample_labels <- fromJSON(sample_labels_path)
head(sample_labels)
class(sample_labels)

# for (row in sample_labels)


print("Breaaaaaak")


train_x = bind_rows(training_data)
# print(training_data)
# print(train_x)

train_y <- logical(nrow(train_x))

for(i in 1:nrow(train_x)){
    row <- train_x[i, ]
    # print(row)
    print(row[["sample_id"]])
    print(sample_labels[[row[["sample_id"]]]])
    print(row[["seg_id"]])
    test <- row[["seg_id"]] %in% sample_labels[[row[["sample_id"]]]]
    print(test)
    print(class(test))
    train_y[i] <- test
}

print(train_y)

model <- train_iimi(train_x, train_y)

predict_iimi(training_data[[1]], "xgb", model)


# test_function <- function (row){
#     sample_labels %>% 
#     filter(virus_id == row[[1]]) %>%
#     nrow() > 0
# }



# train_y <- train_x %>%
#             rowwise() %>%
#             summarise(present =function (row){
#                 sample_labels %>% 
#                 filter(virus_id == row[[1]]) %>%
#                 nrow() > 0
#             })


# test_bool

# df_data = lapply(rle_data, function()))




# for(row in rle_data){
#     row
#     seq_id = row$
# }
 



# bam_file_paths <- lapply(bam_file_paths, function(file_name){
#     file.path("./", bam_path, file_name)
# })



# rle_data <- lapply( bam_file_paths, convert_bam_to_rle)





