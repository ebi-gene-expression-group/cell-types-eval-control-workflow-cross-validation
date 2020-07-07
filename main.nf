#!/usr/bin/env nextflow 

// Workflow to perform k-fold cross validation for evaluation of scRNA-Seq cell predictors over FTP data 

// download data
if(params.download_data.run == "True"){
	
	// extract matrix types required by the tools, conditioned on them being 'on' 
	tool_switch = ["True":0, "False":1]
	garnett_matrix_type = [params.garnett.matrix_type, null][tool_switch[params.garnett.run]]
	scmap_cluster_matrix_type = [params.scmap_cluster.matrix_type, null][tool_switch[params.scmap_cluster.run]]
	scmap_cell_matrix_type = [params.scmap_cell.matrix_type, null][tool_switch[params.scmap_cell.run]]
	scpred_matrix_type = [params.scpred.matrix_type, null][tool_switch[params.scpred.run]]
	
	UNIQUE_MATRIX_TYPES = Channel.of(
		garnett_matrix_type,
		scmap_cluster_matrix_type,
		scmap_cell_matrix_type,
		scpred_matrix_type)
		.filter{ it != null }
		.unique()
	
	// parse txt file with dataset accessions into a queue channel; build combinations with matrix types 
	IMPORT_PARAMS = Channel
                .fromPath(params.download_data.import_datasets)
                .splitCsv(header:false, sep:",")
                .combine(UNIQUE_MATRIX_TYPES)
    
process fetch_query_data{
	publishDir "${params.data_dir}${dataset_id}", mode: 'copy'
	conda "${baseDir}/envs/load_data.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
	maxRetries 5
	memory { 16.GB * task.attempt }
	input:
	tuple val(dataset_id), val(seq_method), val(num_clust), val(barcode_col), val(cell_type_col), val(matrix_type) from IMPORT_PARAMS
	
	output:
	tuple file(dataset_id), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type), val(num_clust) into FETCH_DATA
	"""
	if [ ${seq_method} ==  "droplet" ]; then
        MATRIX_TYPE_UPD="CPM"
    	else
        MATRIX_TYPE_UPD=${matrix_type}
    	fi
    	get_experiment_data.R\
                --accesssion-code ${dataset_id}\
                --config-file ${params.download_data.scxa_import_config_file}\
             	--matrix-type ${matrix_type}\
                --output-dir-name ${dataset_id}\
                --get-sdrf ${params.download_data.get_sdrf}\
                --get-condensed-sdrf ${params.download_data.get_cond_sdrf}\
                --get-idf ${params.download_data.get_idf}\
                --get-marker-genes ${params.download_data.get_marker_genes}\
                --number-of-clusters ${num_clust}
	"""
    }
}

// condensed sdrf files need 'un-melting' 
process unmelt_sdrf_query {
	conda "${baseDir}/envs/exp_metadata.yaml"
	
	input:
	tuple file(data), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type), val(num_clust) from FETCH_DATA
	output:
        tuple file("${data}.${matrix_type}"), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type), val(num_clust) into FETCH_DATA_UNMELT
	file("${data}.${matrix_type}/10x_data/barcodes.tsv") into BARCODES
	val(dataset_id) into DATA_ID
	"""
	unmelt_condensed.R\
		-i ${data}/condensed-sdrf.tsv\
		-o ${data}/unmelted_sdrf.tsv\
		--retain-types ${params.unmelt_sdrf.retain_types}\
		--has-ontology                 
	# rename data dir name to avoid downstream file name collision
	mv ${data} ${data}.${matrix_type} 
	"""
	}

// convert Barcodes to value channel 
BARCODES_FOLDS = BARCODES.first()  
// convert dataset ID to value channel
DATASET_ID = DATA_ID.first()

// generate folds of cell indices for cross validation
process generate_folds{
	publishDir "${params.output_dir}/${dataset_id}", mode: 'copy'
	conda "${baseDir}/envs/r-caret.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	input:
	val(dataset_id) from DATASET_ID
        file(barcodes) from BARCODES_FOLDS

	output:
	file("folds_indices/*") into K_FOLD_CELL_INDEXES 
	"""
	caret-create-folds.R\
	        --input-barcodes ${barcodes}\
	        --k-folds-number ${params.generate_folds.folds_k}\
	        --dataset-id ${dataset_id}\
		--output-dir folds_indices
	"""
}

// flatten channel and set names to 'Fold_X'
K_FOLD_CELL_INDEXES
	.flatten()
	.map{ f -> tuple("${f}".split("\\.")[1], f) }
	.set{FOLDS}

// add fold name to data channels
FOLD_DATA = FOLDS.combine(FETCH_DATA_UNMELT)

// split data into train and test sets based on indices
process split_train_test{
	publishDir "${params.output_dir}/${dataset_id}/split_data", mode: 'copy'
	conda "${baseDir}/envs/r-caret.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	input:
	tuple val(fold), file(indx), file(data), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type), val(num_clust) from FOLD_DATA
	
	output:
	set val(fold), val(dataset_id),  file("${dataset_id}.test.${matrix_type}.${fold}.zip"), file("${dataset_id}.train.${matrix_type}.${fold}.zip") into SPLIT_DATA
	set val(fold), val(dataset_id) into FOLD_N_SPLIT_DATA
	"""
	split-train-test-data.R\
		--input-matrix "${data}/10x_data/matrix.mtx"\
		--input-matrix-type ${matrix_type}\
		--dataset-id ${dataset_id}\
		--input-barcodes "${data}/10x_data/barcodes.tsv"\
		--input-cell-indexes ${indx}\
		--input-unmelt-sdrf "${data}/unmelted_sdrf.tsv"\
		--fold-n ${fold}
	# mv genes and markers files
	parallel cp -v "${data}/10x_data/genes.tsv" ::: ${dataset_id}.test.${matrix_type}.${fold}/10x_data ${dataset_id}.train.${matrix_type}.${fold}/10x_data
	parallel cp -v "${data}/marker_genes_${num_clust}.tsv" ::: ${dataset_id}.test.${matrix_type}.${fold}/ ${dataset_id}.train.${matrix_type}.${fold}
	# zip files by test or train
	zip -r ${dataset_id}.test.${matrix_type}.${fold}.zip ${dataset_id}.test.${matrix_type}.${fold} 
	zip -r ${dataset_id}.train.${matrix_type}.${fold}.zip ${dataset_id}.train.${matrix_type}.${fold} 
	"""
}

// group data by fold (and dataset_id)
GROUPED_DATA = SPLIT_DATA.groupTuple(by:[0, 1])

// move same fold data into single dir
process group_fold_data {
	input:
	tuple val(fold), val(dataset_id), file(x), file(x) from GROUPED_DATA
	output: 
	tuple val(fold), val(dataset_id), file('fold_dir') into RUN_PREDICTORS 
	"""
	mkdir -p fold_dir/
      	for file in '*.zip' 
      	do
              mv \${file} fold_dir
      	done
	"""	

}

// Run cell-types-eval-control-workflow which runs the predictor tools and label analysis
process run_predictors{
	publishDir "${params.output_dir}/${dataset_id}/label_analysis", mode: 'copy'
	conda "${baseDir}/envs/nextflow.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  || task.attempt < 3  ? 'retry' : 'finish' }   
	maxRetries 3
	maxForks 1 
	memory { 16.GB * task.attempt }
	
	input:
	tuple val(fold), val(dataset_id), file(fold_dir) from RUN_PREDICTORS 
	
	output:
	tuple file("${dataset_id}.${fold}.${params.run_predictors.tool_perf_table}"), file("${dataset_id}.${fold}.${params.run_predictors.tool_perf_pvals}") into LABEL_ANALYSIS
	val(dataset_id) into DATASET_ID_A
	"""
	RESULTS_DIR="\$PWD"
	# launch feature branch from control workflow	
	pushd $CONTROL_WORKFLOW > /dev/null
	git checkout feature/diff_matrix_download > /dev/null
	popd > /dev/null

	nextflow run $CONTROL_WORKFLOW/main.nf\
		-profile ${params.profile}\
		-c "${baseDir}/nextflow.config"\
		--label_analysis_outdir \$RESULTS_DIR\
		--input_data "${fold_dir}"
	# rename prediction outputs to include dataset ID and fold value
	mv ${params.run_predictors.tool_perf_table} ${dataset_id}.${fold}.${params.run_predictors.tool_perf_table}
	mv ${params.run_predictors.tool_perf_pvals} ${dataset_id}.${fold}.${params.run_predictors.tool_perf_pvals}
	"""
} 

// combine results into a single directory
process combine_results{
	input:
	file(folds_label) from LABEL_ANALYSIS.collect()
	
	output:
	file('results_dir') into COMBINED_RESULTS
	
	"""
	mkdir -p results_dir/
	for file in ${folds_label}
	do
		mv \${file} results_dir
	done
	"""
}

// average folds output
process avg_folds {
	publishDir "${params.output_dir}/${dataset_id}/final_output", mode: 'copy'
	conda "${baseDir}/envs/r-caret.yaml"
	errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
	maxRetries 5
	memory { 16.GB * task.attempt }
	
	input:
	file ('results_dir') from COMBINED_RESULTS
	val(dataset_id) from DATASET_ID_A 
	output:
	    file("${dataset_id}.avg.tool_perf_table.tsv") into AVG_TOOL_PERF_TABLE 
	    file("${dataset_id}.avg.tool_perf_pvals.tsv") into AVG_TOOL_PERF_PVALS
	"""
	avg-tables-cell-label-analysis.R\
		--input-dir ${results_dir}\
	        --dataset-id ${dataset_id}\
		--tool-perf-table ${params.run_predictors.tool_perf_table}\
		--tool-perf-pvals ${params.run_predictors.tool_perf_pvals}
	"""
}
