#!/usr/bin/env Rscript 

#Average performance tables of the label analysis corresponding to the different folds into a single one.
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-t", "--input-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Diectory path to output of previous processes where files are stored by Fold.'
  ),
  make_option(
    c("-d", "--dataset-id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Dataset ID of the dataset used'
  ), 
  make_option(
    c("-a", "--avg-tool-perf-table"),
    action = "store",
    default = "avg_tool_perf_table.tsv",
    type = 'character',
    help = 'Cross-fold averaged tool performance table.'
  ),
  make_option(
    c("-v", "--avg-tool-perf-pvals"),
    action = "store",
    default = "avg_tool_perf_pvals.tsv",
    type = 'character',
    help = 'Cross-fold averaged tool p values table.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_dir"))

if(!file.exists(opt$input_dir)) stop("Input directory containing label-analysis outputs does not exist.")

#read files
tool_perf_tables_paths <- list.files(opt$input_dir, pattern="tool_perf_table", full.names=T)
tool_perf_table_list <- lapply(tool_perf_tables_paths, FUN = function(x) read.table(x, header = T))

tool_perf_tables_paths <- list.files(opt$input_dir, pattern="tool_perf_pvals", full.names=T)
tool_perf_pvals_list <- lapply(tool_perf_tables_paths, FUN = function(x) read.table(x, header = T))

#merge tables
merge_tool_perf_tables <- Reduce(rbind, tool_perf_table_list)
merge_tool_perf_pvals <- Reduce(rbind, tool_perf_pvals_list)
#calculate mean values by tool
mean_tool_perf_tables <- aggregate(x = subset(merge_tool_perf_tables, select = -c(Tool)), by = list(merge_tool_perf_tables$Tool), FUN = mean)
mean_tool_perf_pvals <- aggregate(x = subset(merge_tool_perf_pvals, select = -c(Tool)), by = list(merge_tool_perf_pvals$Tool), FUN = mean)
#replace Group.1 colname to Tool 
colnames(mean_tool_perf_tables) <- replace(colnames(mean_tool_perf_tables), c(1), c("Tool"))
colnames(mean_tool_perf_pvals) <- replace(colnames(mean_tool_perf_pvals), c(1), c("Tool"))
#save tables

write.table(mean_tool_perf_tables, file=paste0(opt$dataset_id, "_", opt$avg_tool_perf_table), sep = "\t")
write.table(mean_tool_perf_pvals, file=paste0(opt$dataset_id, "_", opt$avg_tool_perf_pvals), sep = "\t")
