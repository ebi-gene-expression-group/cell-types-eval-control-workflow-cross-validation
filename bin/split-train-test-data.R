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
    help = 'Path to the input cell indexes to perform test/train data subseting.'
  ),
  make_option(
    c("-b", "--input-barcodes-tsv"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the barcodes tsv to subset based on input_cell_indexes.'
  ),
  make_option(
    c("-f", "--input-genes-tsv"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the genes tsv.'
  ),
  make_option(
    c("-s", "--input-processed-sdrf"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the processed SDRF file.'
  ),
 make_option(
    c("-q", "--output-test-dir"),
    action = "store",
    default = 'test',
    type = 'character',
    help = "Name the test matrix split." 
  ),
 make_option(
    c("-t", "--output-train-dir"),
    action = "store",
    default = 'train',
    type = 'character',
    help = "Name of the train matrix split." 
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_matrix", "input_cell_indexes","input_barcodes_tsv", "input_genes_tsv", "input_processed_sdrf"))

# read 10X sparse matrix
suppressPackageStartupMessages(require(Matrix))
if(!file.exists(opt$input_matrix)) stop("Input 10X matrix file not provided.")
input_matrix <- readMM(opt$input_matrix)
# read input cell indexes 
if(!file.exists(opt$input_cell_indexes)) stop("Input cell labels file does not exist.")
test_cell_indexes <- readRDS(opt$input_cell_indexes)
# read input barcodes 
if(!file.exists(opt$input_barcodes_tsv)) stop("Input barcodes file does not exist.")
barcodes <- read.table(opt$input_barcodes_tsv, header = F, sep = "\t", quote = "")
# read input genes 
if(!file.exists(opt$input_genes_tsv)) stop("Input genes file does not exist.")
genes <- read.table(opt$input_genes_tsv, header = F, sep = "\t", quote = "")
# read input processd srdf 
if(!file.exists(opt$input_processed_sdrf)) stop("Input sdrf file does not exist.")
sdrf <- read.table(opt$input_processed_sdrf, header = T, sep = "\t", quote = "")

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

# Make output directories
dir.create('test')
dir.create('train')

# save matrices
writeMM(test_matrix, file = file.path(opt$output_test_dir,"matrix.mtx"))
writeMM(train_matrix, file = file.path(opt$output_train_dir,"matrix.mtx"))
# save barcodes
write.table(test_barcodes, file = file.path(opt$output_test_dir,"barcodes.tsv"), sep = "\t", quote = F, row.names=F , col.names = F)
write.table(train_barcodes, file = file.path(opt$output_train_dir,"barcodes.tsv"), sep = "\t", quote = F, row.names=F , col.names = F)
# save genes
write.table(genes, file = file.path(opt$output_test_dir,"genes.tsv"), sep = "\t", quote = F, row.names=F , col.names = F)
write.table(genes, file = file.path(opt$output_train_dir,"genes.tsv"), sep = "\t", quote = F, row.names=F , col.names = F)
# save metadata 
write.table(test_sdrf, file = file.path(opt$output_test_dir,"test_sdrf.tsv"), sep = "\t", quote = F,  row.names=F )
write.table(train_sdrf, file = file.path(opt$output_train_dir,"train_sdrf.tsv"), sep = "\t", quote = F,  row.names=F )
