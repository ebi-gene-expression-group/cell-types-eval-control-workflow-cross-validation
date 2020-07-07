#!/usr/bin/env Rscript 

#Generate k-fold cell indexes of input SCE object. 
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-barcodes"),
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
    c("-d", "--dataset-id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Dataset ID of the dataset used'
  ), 
  make_option(
    c("-o", "--output-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the each of the cell indexes for each fold in rds format.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_barcodes", "output_dir"))

# read SCE object
if(!file.exists(opt$input_barcodes)) stop("Input file does not exist.")
barcodes <- read.table(file=opt$input_barcodes, sep = "\t", header=T)
# create Folds with cell indexes
suppressPackageStartupMessages(require(caret))
# list of folds with cell indexes
barcodes_index_list <- createFolds(y = 1:nrow(barcodes), k = opt$k_folds_number)
print(barcodes_index_list)
# save fold's barcode indices
dir.create(opt$output_dir)
for (i in 1:length(barcodes_index_list)){saveRDS(barcodes_index_list[[i]], file=file.path(opt$output_dir, paste0(opt$dataset_id, ".fold_", i, ".indices.rds")))}
