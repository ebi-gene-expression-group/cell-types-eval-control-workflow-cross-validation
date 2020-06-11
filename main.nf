#!/usr/bin/env nextflow 

// Workflow to perform k-fold cross validation over a downloaded dataset (via accession id) or an input dataset

// download query data
if(params.download_data.run == "True"){
    query_data = params.download_data.query_output_dir
    query_n_clust = params.download_data.query_num_clust.toString()
    query_markers = "marker_genes_" + query_n_clust + ".tsv"

process fetch_query_data{
	publishDir "${baseDir}/data", mode: 'copy'
	conda "${baseDir}/envs/load_data.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	output:
	file("${query_data}/query_10x_data") into QUERY_10X_DIR
	file("${query_data}/query_sdrf.txt") into CONDENSED_SDRF
	file("${query_data}/query_${query_markers}") into MARKERS
	"""
	get_experiment_data.R\
		--accesssion-code ${params.download_data.query_acc_code}\
		--expr-data-type ${params.download_data.expr_data_type}\
		--normalisation-method ${params.download_data.normalisation_method}\
		--output-dir-name ${params.download_data.query_output_dir}\
		--get-condensed-sdrf\
		--get-marker-genes\
		--number-of-clusters ${params.download_data.query_num_clust}
	# rename files to avoid name collisions in subsequent processes
	mv ${query_data}/10x_data ${query_data}/query_10x_data
	mv ${query_data}/condensed-sdrf.tsv ${query_data}/query_sdrf.txt
	mv ${query_data}/${query_markers} ${query_data}/query_${query_markers}
	"""
    }
    
}else{
	QUERY_10X_DIR = Channel.fromPath(params.data_10X, checkIfExists: true).first()
	CONDENSED_SDRF = Channel.fromPath(params.condensed_sdrf, checkIfExists: true).first()
	MARKERS = Channel.fromPath(params.markers_metadata, checkIfExists: true).first()
}	

// condensed sdrf files need 'un-melting' 
process unmelt_sdrf_query {
	conda "${baseDir}/envs/exp_metadata.yaml"
	
	input:
	file(condensed_sdrf) from CONDENSED_SDRF
	output:
	file("query_metadata.tsv") into UNMELTED_SDRF
	"""
	unmelt_condensed.R\
		-i ${condensed_sdrf}\
		-o query_metadata.tsv\
		--retain-types ${params.unmelt_sdrf.retain_types}\
		--has-ontology                 
	"""
	}

// generate folds of cell indices for cross validation
process generate_folds{
	publishDir "${params.output_dir}", mode: 'copy'
	conda "${baseDir}/envs/r-caret.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	input:
	file(data_dir) from QUERY_10X_DIR 
	output:
	file("folds_indices/*") into K_FOLD_CELL_INDEXES 
	"""
	caret-create-folds.R\
	        --input-barcodes-tsv "${data_dir}/barcodes.tsv"\
	        --k-folds-number ${params.generate_folds.folds_k}\
	        --output-dir folds_indices
	"""
}

K_FOLD_CELL_INDEXES
    .flatten()
    .map{ f -> tuple("${f.simpleName}", f) }
    .set{FOLDS_BY_K}
	
// split data into train and test sets based on indices
process split_train_test{
	//publishDir "${params.output_dir}/split_data", mode: 'copy'
	conda "${baseDir}/envs/r-caret.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	input:
	    set val(fold), file(k_fold_cell_indexes) from FOLDS_BY_K 
	    file(data_dir) from QUERY_10X_DIR 
	    file(sdrf) from UNMELTED_SDRF
	output:
	    set val("${fold}"), file("test.zip"), file("train.zip") into SPLIT_DATA 
	"""
	split-train-test-data.R\
		--input-matrix "${data_dir}/matrix.mtx"\
		--input-barcodes-tsv "${data_dir}/barcodes.tsv"\
		--input-genes-tsv "${data_dir}/genes.tsv"\
		--input-cell-indexes ${k_fold_cell_indexes}\
		--input-processed-sdrf ${sdrf}\
		--output-test-dir test\
		--output-train-dir train
	#zip files by test or train
	zip -r test.zip test
	zip -r train.zip train
	"""
}

// Run cell-types-eval-control-workflow which runs the predictor tools
process run_cell_types_eval {
	publishDir "${params.output_dir}/label_analysis", mode: 'copy'
	conda "${baseDir}/envs/nextflow.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  || task.attempt < 3  ? 'retry' : 'finish' }   
	maxRetries 3
	maxForks 1 
	memory { 16.GB * task.attempt }
	
	input:
	    tuple val(fold), file(test), file(train) from SPLIT_DATA
	    file(marker_genes) from MARKERS
	output:
	    set val(fold), file("tool_perf_table.tsv") into TOOL_PERF_TABLE
	    set val(fold), file("tool_perf_pvals.tsv") into TOOL_TABLE_PVALS
	"""
	test_zipdir=\$(unzip -qql $test | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/.*||')
	unzip $test
	test_sdrf=\${test_zipdir}/test_sdrf.tsv\
	train_zipdir=\$(unzip -qql $train | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/.*||')
	unzip $train
	train_sdrf=\${train_zipdir}/train_sdrf.tsv\
	RESULTS_DIR="\$PWD"
  
	# launch develop branch from control workflow	
	pushd $CONTROL_WORKFLOW > /dev/null
	git checkout develop > /dev/null
	popd > /dev/null

	nextflow run $CONTROL_WORKFLOW/main.nf\
		-profile ${params.profile}\
		-c "${baseDir}/nextflow.config"\
		--label_analysis_outdir \$RESULTS_DIR\
		--ref_10x_dir \${train_zipdir}\
		--query_10x_dir \${test_zipdir}\
		--unmelt_sdrf_ref \${train_sdrf}\
		--unmelt_sdrf_query \${test_sdrf}\
		--ref_markers ${marker_genes}
	"""
} 

// join performance and p-value tables into a single channel
OUT_CH = TOOL_PERF_TABLE.join(TOOL_TABLE_PVALS)

// rename files to include fold number
process rename_results{
	publishDir "${params.output_dir}/${fold}", mode: 'copy'
  
	input:
	set val(fold), file(tool_perf_table), file(tool_perf_pvals) from OUT_CH
	output:
	set file("${fold}_${tool_perf_table}"), file("${fold}_${tool_perf_pvals}") into RENAMED_RESULTS
	
	"""
	mv ${tool_perf_table} ${fold}_${tool_perf_table}
	mv ${tool_perf_pvals} ${fold}_${tool_perf_pvals}
	"""
}

// combine results into a single directory
process combine_results{
	publishDir "${params.output_dir}", mode: 'copy'
	
	input:
	file(folds_output) from RENAMED_RESULTS.collect()
	output:
	file('results_dir') into COMBINED_RESULTS_DIR
	"""
	mkdir -p results_dir/
	for file in ${folds_output}
	do
		mv \${file} results_dir
	done
	"""
}

// average folds output
process average_folds_output {
	publishDir "${params.output_dir}/final_output", mode: 'copy'
	conda "${baseDir}/envs/r-caret.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	input:
	file ('results_dir') from COMBINED_RESULTS_DIR
	output:
	    file("avg_tool_perf_table.tsv") into AVG_TOOL_PERF_TABLE 
	    file("avg_tool_perf_pvals.tsv") into AVG_TOOL_PERF_PVALS
	"""
	avg-tables-cell-label-analysis.R\
		--input-dir ${results_dir}\
		--avg-tool-perf-table ${params.average_folds_output.avg_tool_perf_table}\
		--avg-tool-perf-pvals ${params.average_folds_output.avg_tool_perf_pvals}
	"""
}
