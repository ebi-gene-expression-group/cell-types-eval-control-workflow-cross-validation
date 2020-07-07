#!/usr/bin/env Rscript 

#This script is to perform subsetting of dataset based on the cel labels generated in the k-fold cv previous process. 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
 make_option(
    c("-i", "--input-matrix"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to the 10X matrix to split into test and train subsets."
  ), 
 make_option(
    c("-y", "--input-matrix-type"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Type of input matrix e.g 'CPM', 'TPM', 'RAW'"
  ), 
 make_option(
    c("-d", "--dataset-id"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Dataset ID."
  ), 
  make_option(
    c("-c", "--input-cell-indexes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input cell indexes to perform test/train data subseting.'
  ),
  make_option(
    c("-b", "--input-barcodes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the barcodes tsv to subset based on input_cell_indexes.'
  ),
  make_option(
    c("-s", "--input-unmelt-sdrf"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the unmelt SDRF file.'
  ),
  make_option(
    c("-f", "--fold-n"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Fold number' 
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_matrix", "input_matrix_type","dataset_id", "input_cell_indexes","input_barcodes", "input_unmelt_sdrf", "fold_n"))

# read input matrix 
suppressPackageStartupMessages(require(Matrix))
if(!file.exists(opt$input_matrix)) stop("Input 10X CPM matrix file not provided.")
input_matrix <- readMM(opt$input_matrix)
# args
matrix_type <- opt$input_matrix_type
dataset_id <- opt$dataset_id
fold_n <- opt$fold_n
# read input cell indexes 
if(!file.exists(opt$input_cell_indexes)) stop("Input cell labels file does not exist.")
test_cell_indexes <- readRDS(opt$input_cell_indexes)
# read input barcodes 
if(!file.exists(opt$input_barcodes)) stop("Input barcodes file does not exist.")
barcodes <- read.table(opt$input_barcodes, header = F, sep = "\t", quote = "")
# read input processd srdf 
if(!file.exists(opt$input_unmelt_sdrf)) stop("Input sdrf file does not exist.")
sdrf <- read.table(opt$input_unmelt_sdrf, header = T, sep = "\t", quote = "")
# order sdrf file with barcodes.tsv file order
sdrf <- sdrf[order(match(sdrf[,"id"], barcodes[,1]), decreasing = F), ]
# generate train indexes
train_cell_indexes <- c(1:ncol(input_matrix))[!(1:ncol(input_matrix) %in% test_cell_indexes)]

# split data into test and train sets
test_matrix <- input_matrix[, test_cell_indexes]
test_barcodes <- data.frame(barcodes[test_cell_indexes, ])
test_sdrf <- data.frame(sdrf[ test_cell_indexes, ])

train_matrix <- input_matrix[, train_cell_indexes]
train_barcodes <- data.frame(barcodes[train_cell_indexes, ]) 
train_sdrf <- data.frame(sdrf[train_cell_indexes, ])

# compose file name 
test_dir <- paste0(dataset_id, ".test.", matrix_type, ".", fold_n)
train_dir <- paste0(dataset_id, ".train.", matrix_type, ".", fold_n)
dir.create(test_dir)
dir.create(paste0(test_dir, "/10x_data"))
dir.create(train_dir)
dir.create(paste0(train_dir, "/10x_data"))

# save matrices
writeMM(test_matrix, file = file.path(test_dir,"10x_data", "matrix.mtx"))
writeMM(train_matrix, file = file.path(train_dir, "10x_data", "matrix.mtx"))
# save barcodes
write.table(test_barcodes, file = file.path(test_dir, "10x_data", "barcodes.tsv"), sep = "\t", quote = F, row.names=F , col.names = F)
write.table(train_barcodes, file = file.path(train_dir, "10x_data", "barcodes.tsv"), sep = "\t", quote = F, row.names=F , col.names = F)
# save metadata 
write.table(test_sdrf, file = file.path(test_dir, "unmelted_sdrf.tsv"), sep = "\t", quote = F,  row.names=F )
write.table(train_sdrf, file = file.path(train_dir, "unmelted_sdrf.tsv"), sep = "\t", quote = F,  row.names=F )
