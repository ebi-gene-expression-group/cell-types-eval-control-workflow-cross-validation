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
    c("-c", "--input-cell-indexes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input cell indexes to perform test/train data subseting'
  ),
  make_option(
    c("-b", "--input-barcodes-tsv"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the barcodes tsv to subset based on input_cell_indexes.'
  ),
  make_option(
    c("-f", "--input-features-tsv"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the features tsv.'
  ),
 make_option(
    c("-q", "--output-test-dir"),
    action = "store",
    default = 'test',
    type = 'character',
    help = "Name the test matrix split." #where cell identity will be inferred based on the model generated over train sce split."
  ),
 make_option(
    c("-t", "--output-train-dir"),
    action = "store",
    default = 'train',
    type = 'character',
    help = "Name of the train matrix split." #where the the model for inferring cell identity will be generated."
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_matrix", "input_cell_indexes","input_barcodes_tsv", "input_features_tsv"))

#read 10X sparse matrix
suppressPackageStartupMessages(require(Matrix))
if(!file.exists(opt$input_matrix)) stop("Input 10X matrix file not provided.")
input_matrix <- readMM(opt$input_matrix)
#read input cell indexes 
if(!file.exists(opt$input_cell_indexes)) stop("Input cell labels file does not exist.")
test_cell_indexes <- readRDS(opt$input_cell_indexes)
#read input barcodes 
if(!file.exists(opt$input_barcodes_tsv)) stop("Input barcodes file does not exist.")
barcodes <- read.table(opt$input_barcodes_tsv, header = F, sep = "\t", quote = "")
#read input features 
if(!file.exists(opt$input_features_tsv)) stop("Input features file does not exist.")
features <- read.table(opt$input_features_tsv, header = F, sep = "\t", quote = "")

train_cell_indexes <- c(1:ncol(input_matrix))[!(1:ncol(input_matrix) %in% test_cell_indexes)]

#split data into test and train
test_matrix <- input_matrix[, test_cell_indexes]
train_matrix <- input_matrix[, train_cell_indexes]

test_barcodes <- data.frame(barcodes[test_cell_indexes, ])
train_barcodes <- data.frame(barcodes[train_cell_indexes, ]) 

#save gzipped data


#save data individually
writeMM(test_matrix, file = paste0(opt$output_test_dir,"matrix.mtx"))
writeMM(train_matrix, file = paste0(opt$output_train_dir,"matrix.mtx"))

write.table(test_barcodes, file = paste0(opt$output_test_dir,"barcodes.tsv"), sep = "\t", quote = F, col.names = F)
write.table(train_barcodes, file = paste0(opt$output_train_dir,"barcodes.tsv"), sep = "\t", quote = F, col.names = F)

write.table(features, file = paste0(opt$output_test_dir,"features.tsv"), sep = "\t", quote = F, col.names = F)
write.table(features, file = paste0(opt$output_train_dir,"features.tsv"), sep = "\t", quote = F, col.names = F)

#test_mat_list <- lapply(test_cell_indexes, function(fold) mat[, fold])
#test_barcodes_list <- lapply(test_cell_indexes, function(fold) data.frame(barcodes[fold,]))
#train_mat_list <- lapply(test_cell_indexes, function(fold) y$mat[, !(1:ncol(mat) %in% fold)])
#train_barcodes_list <- lapply(test_cell_indexes, function(fold) data.frame(barcodes[!(1:ncol(mat) %in% fold), ]))
#
##save TEST Matrices and Barcodes
#mapply(function(X, Y){writeMM(X, file=paste0(opt$output_test_dir, ".", Y, ".matrix.mtx"))}, X= test_mat_list, Y=as.list(names(test_mat_list)))
#mapply(function(X, Y){write.table(X, file=paste0(opt$output_test_dir, ".", Y, ".barcodes.tsv"), sep = "\t", quote = F, col.names = F)}, X= test_barcodes_list, Y=as.list(names(test_barcodes_list)))
#
##save TRAIN Matrices and Barcodes
#mapply(function(X, Y){writeMM(X, file=paste0(opt$output_train_dir, ".", Y, ".matrix.mtx"))}, X= train_mat_list, Y=as.list(names(train_mat_list)))
#mapply(function(X, Y){write.table(X, file=paste0(opt$output_train_dir, ".", Y, ".barcodes.tsv"), sep = "\t", quote = F, col.names = F)}, X= train_barcodes_list, Y=as.list(names(train_barcodes_list)))

#split data into test and train
#test_mat_list <- lapply(test_cell_indexes, function(fold) mat[, fold])
#test_barcodes_list <- lapply(test_cell_indexes, function(fold) data.frame(barcodes[fold,]))
#train_mat_list <- lapply(test_cell_indexes, function(fold) mat[, !(1:ncol(mat) %in% fold)])
#train_barcodes_list <- lapply(test_cell_indexes, function(fold) data.frame(barcodes[!(1:ncol(mat) %in% fold), ]))
#
##save TEST Matrices and Barcodes
#mapply(function(X, Y){writeMM(X, file=paste0(opt$test_name, ".", Y, ".matrix.mtx"))}, X= test_mat_list, Y=as.list(names(test_mat_list)))
#mapply(function(X, Y){write.table(X, file=paste0(opt$test_name, ".", Y, ".barcodes.tsv"), sep = "\t", quote = F, col.names = F)}, X= test_barcodes_list, Y=as.list(names(test_barcodes_list)))
#
##save TRAIN Matrices and Barcodes
#mapply(function(X, Y){writeMM(X, file=paste0(opt$train_name, ".", Y, ".matrix.mtx"))}, X= train_mat_list, Y=as.list(names(train_mat_list)))
#mapply(function(X, Y){write.table(X, file=paste0(opt$train_name, ".", Y, ".barcodes.tsv"), sep = "\t", quote = F, col.names = F)}, X= train_barcodes_list, Y=as.list(names(train_barcodes_list)))