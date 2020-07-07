# cell-types-eval-control-workflow-cross-validation


A k-fold cross validation workflow for evaluating the performance of scRNA-Seq cell annotation tools. 

The input dataset is split into k groups, each of which is treated as test data set, and the remaining as train data set in the nested 
[cell-types-eval-control-workflow](https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow-cross-validation/tree/feature/integrate_celltypes_workflow) which contains the following cell type predictors: 
* [garnett-eval-workflow](https://github.com/ebi-gene-expression-group/garnett-eval-workflow)
* [scmap-eval-workflow](https://github.com/ebi-gene-expression-group/scmap-eval-workflow) - 2 versions: 'cluster' and 'cell'
* [scpred-eval-workflow](https://github.com/ebi-gene-expression-group/scpred-eval-workflow)
* [label-analysis-eval-workflow](https://github.com/ebi-gene-expression-group/label-analysis-eval-workflow)

Workflow diagram:

![](https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow-cross-validation/blob/develop/cross-validation-pipeline.png)

### Pipeline configuration

##### Number of folds (k) 
The number of divisions in which to split the input data can be configured via the `generate_folds.folds_k` parameter in the config file.

##### Enable/disable cell predictors
Cell prediction tools can be custom enabled and disabled from the `nextflow.config` file. Additionally, matrix type required by each tool is specified here, allowing data import step to fetch different data matrices (TPMs, CPMs, raw or filtered) to be fetched from the server. 

##### Import datasets
Data to be imported from the FTP server is specified in the config file via a [CSV file](https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow-cross-validation/blob/develop/datasets.txt) containing the following fields: 
* `dataset id` (from SCXA)
* `technology type` _("droplet" or "smart-seq")_
* `matrix type` _(raw, filtered, CPM- or TPM-normalised)_
* `number of clusters in marker gene file` 
* `barcode column` (in SDRF file) 
* `cell type column` 

### Triggering the pipeline
You will need Nextflow installed to run analyses. The best way to install it is via Conda. 
To initially create a clean environment and avoid potential dependency conflicts run the following commands:

```
conda create --name nextflow
conda activate nextflow
conda install nextflow
```

To pull the latest version of the cross-validation workflow and it's nested components issue the following command from the `cell-types-eval-control-workflow-cross-validation` directory:
```
./bin/fetch-control-and-tool-eval-workflows.sh
```
Finally, to run the pipeline run: 
```
nextflow run main.nf --profile cluster
```
In case execution is halted, you can always resume the execution from the point it stopped, specifying: 
```
nextflow run main.nf --profile cluster -resume
```

### Results

The final output of the workflow are 2 tables containing the metrics of the label analysis of the cell predictions averaged for all folds (as defined by `fold_k`), and it's associated significances. See an final output example [here](https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow-cross-validation/tree/master/example_output).

These, and additional outputs (fold data splitting, individual fold label analysis) are stored in the directory `output_dir/dataset_id` being `output_dir` defined in the config file. 


### Example Run on dataset E-ENAD-27
Here is an example on how to run the pipeline on SCXA dataset `E-ENAD-27`. 

First of all we must specify the `dataset id`, `technology type`, `matrix type`, `number of clusters in marker gene file`, `barcode column`, and `cell type column` in a CSV file (located in the same directory where the pipeline is going to be run): 
```
E-ENAD-27,smart-seq,7,id,inferred.cell.type
```

Next, running the following command will **pull the latest version** of the cross-validation workflow and it's nested components:
```
./bin/fetch-control-and-tool-eval-workflows.sh
```

Once the workflows are updated, we can **edit the `nextflow.config`** to: 
- Specify the output directory: `output_dir`
- Edit the metadata fields of the imported SDRF file.
- Specify wether we want SDRF file to be unmelt: `unmelt_sdrf.run`.
- Specify the number of folds to split the data to perform cross validation: `generate_folds.folds_k`
- Enable/disable the different cell predictors: `scpred.run`, `scmap_cluster.run` ...
- Specify the matrix type required by each cell predictor: `scpred.matrix_type`,  `scmap_cluster.matrix_type`...

**Note:** Current values of matrix type for each tool are those recommended by the authors.
- Specify different parameters of the cell predictors, for instance: `scpred.model == 'svmRadialWeights'`
- Enable/disable label analysis with the cell type predictions done by the methods.

Finally we **run the pipeline** by issuing the following command: 
```
nextflow run main.nf --profile cluster
```

Once the execution is finalized we will find the following **outputs** in a directory named with the `dataset id` within the defined `output_dir`: 
- `/fold_indices`: Containing `.rds` files with the cell indexes for each of the folds.
- `/split_data`: Containing the test-train data split for each fold based. 
- `/label_analysis`: Containing the results of the label analysis for each fold. 
- `/final_output`: Containing the all-fold's average results tables.

