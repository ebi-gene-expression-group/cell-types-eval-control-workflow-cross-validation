#!/usr/bin/env nextflow 

// Workflow to perform k-fold cross validation over the imput dataset

//cross validation
if(params.run == "True"){
	// generate folds of cell indexes
	
K_RANGE = Channel.from(1..params.generate_folds.k_folds_num).map{it -> "Fold_$it"} 
	
	process generate_folds{
		publishDir "${params.output_dir}", mode: 'copy'
        	conda "${baseDir}/envs/r-caret.yaml"

		errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
        	maxRetries 5
        	memory { 16.GB * task.attempt }
		//input:
        	//val k from K_RANGE
		output:
        	//tuple val(k), file("output_cell_indexes.rds") into K_FOLD_CELL_INDEXES 
        	file("*.output_cell_indexes.rds") into K_FOLD_CELL_INDEXES 
        	
		script:
		"""
        	caret-create-folds.R\
        	        --input-barcodes-tsv ${params.barcodes}\
        	        --k-folds-number ${params.generate_folds.k_folds_num}\
        	        --output-cell-indexes ${params.generate_folds.output_cell_indexes}
        	"""
	}

//NOTE: I did this (faltten and map) as I thought the reason for the next process not being iterative was that files were being overwritten, but doing this does not solve anything. Split data process stil runs all in once. WWHY!?
FOLD_INDEXES = K_FOLD_CELL_INDEXES.flatten().map{f -> tuple("${f}".split("\\.")[1], f) }

// Data value channels 
MATRIX_CH = Channel.fromPath(params.matrix)
BARCODES_CH = Channel.fromPath(params.barcodes)
FEATURES_CH = Channel.fromPath(params.features)
SDRF_CH = Channel.fromPath(params.sdrf_processed)

//PROCESS NOT ITERATING!! WHY?
// split data based on cell index folds
	process split_train_test{
	// We're hard coding the output dir in the config file	
		publishDir "${params.output_dir}", mode: 'copy'
		conda "${baseDir}/envs/r-caret.yaml"
		
		errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
                maxRetries 5
                memory { 16.GB * task.attempt }

                input:
		tuple val(fold), file(k_fold_cell_indexes) from FOLD_INDEXES 
		file(matrix) from MATRIX_CH 
		file(barcodes) from BARCODES_CH 
		file(features) from FEATURES_CH 
		file(sdrf) from SDRF_CH
		
		output:
		tuple val(fold), file("test.zip"), file("train.zip") into SPLIT_DATA 
		
		"""
                rm -rf ${params.split_train_test.output_test_dir} && mkdir ${params.split_train_test.output_test_dir}
                rm -rf ${params.split_train_test.output_train_dir} && mkdir ${params.split_train_test.output_train_dir}
                #mkdir ${params.split_train_test.output_train_dir}

		split-train-test-data.R\
                        --input-matrix ${matrix}\
                        --input-cell-indexes ${k_fold_cell_indexes}\
			--input-barcodes-tsv ${barcodes}\
			--input-features-tsv ${features}\
			--input-processed-sdrf ${sdrf}\
                        --output-test-dir ${params.split_train_test.output_test_dir}\
                        --output-train-dir ${params.split_train_test.output_train_dir}
		#zip files by test or train
		zip -r test.zip ${params.split_train_test.output_test_dir}
		zip -r train.zip ${params.split_train_test.output_train_dir}
                """
        }


SPLIT_DATA.subscribe{ println it } 
}
