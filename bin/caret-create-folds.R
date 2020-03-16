#!/usr/bin/env Rscript 

#Generate k-fold cell indexes of input SCE object. 
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-barcodes-tsv"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the barcodes tsv where the folds are going to be computed.'
  ),
  make_option(
    c("-k", "--k-folds-number"),
    action = "store",
    default = 5,
    type = 'integer',
    help = 'Number of groups to split the data in.'
  ), 
  make_option(
    c("-o", "--output-cell-indexes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the each of the cell indexes for each fold in rds format.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_barcodes_tsv", "output_cell_indexes"))

#read SCE object
if(!file.exists(opt$input_barcodes_tsv)) stop("Input file does not exist.")
#sce <- readRDS(opt$input_sce_object)
barcodes <- read.table(file=opt$input_barcodes_tsv, sep = "\t", header=T)
#if object has no colnames, append them
#if(is.null(colnames(sce))) colnames(sce) <- sce$Barcode

#create Folds with cell indexes
suppressPackageStartupMessages(require(caret))
#list of folds with cell indexes
barcodes_index_list <- createFolds(y = 1:nrow(barcodes), k = opt$k_folds_number)
print(barcodes_index_list)

for (i in 1:length(barcodes_index_list)){saveRDS(barcodes_index_list[[i]], file= paste0("index.Fold_", i, ".", opt$output_cell_indexes))}
