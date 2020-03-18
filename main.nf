#!/usr/bin/env nextflow 

// Workflow to perform k-fold cross validation over the imput dataset

K_RANGE = Channel.from(1..params.generate_folds.k_folds_num).map{it -> "Fold_$it"} 
MATRIX_CH = Channel.fromPath(params.matrix).first()
SDRF_CH = Channel.fromPath(params.sdrf_processed).first()
	
process generate_folds{
    publishDir "${params.output_dir}", mode: 'copy'
        conda "${baseDir}/envs/r-caret.yaml"

    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
        maxRetries 5
        memory { 16.GB * task.attempt }
    
    input:
        file(matrix) from MATRIX_CH        
    
    output:
        file("folds/*") into K_FOLD_CELL_INDEXES 

    script:
    """
        zipdir=\$(unzip -qql $matrix | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/.*||')
        unzip -p $matrix \${zipdir}/barcodes.tsv > barcodes.tsv
        
        caret-create-folds.R\
                --input-barcodes-tsv barcodes.tsv\
                --k-folds-number ${params.generate_folds.k_folds_num}\
                --output-dir folds
        """
}

K_FOLD_CELL_INDEXES
    .flatten()
    .map{ f -> tuple("${f.simpleName}", f) }
    .set{FOLDS_BY_K}
	
process split_train_test{
    publishDir "${params.output_dir}", mode: 'copy'
    conda "${baseDir}/envs/r-caret.yaml"
    
    errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }
    maxRetries 5
    memory { 16.GB * task.attempt }

    input:
        set val(fold), file(k_fold_cell_indexes) from FOLDS_BY_K 
        file(matrix) from MATRIX_CH 
        file(sdrf) from SDRF_CH
    
    output:
        set val(fold), file("test.zip"), file("train.zip") into SPLIT_DATA 
    
    """
    zipdir=\$(unzip -qql $matrix | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/.*||')
    unzip $matrix
    
    split-train-test-data.R\
        --input-matrix \${zipdir}/matrix.mtx\
        --input-barcodes-tsv \${zipdir}/barcodes.tsv \
        --input-features-tsv \${zipdir}/genes.tsv\
        --input-cell-indexes ${k_fold_cell_indexes}\
        --input-processed-sdrf ${sdrf}\
        --output-test-dir test\
        --output-train-dir train
    
    #zip files by test or train
    zip -r test.zip test
    zip -r train.zip train
    """
}

SPLIT_DATA.view{ it }
