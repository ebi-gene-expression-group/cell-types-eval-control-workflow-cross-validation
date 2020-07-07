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
