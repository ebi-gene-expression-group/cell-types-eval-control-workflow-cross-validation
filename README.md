# cell-types-eval-control-workflow-cross-validation

A k-fold cross validation workflow for evaluation of cell type classification tools consistency. Input dataset is split into k groups, each of which is treated as test data set, and the remaining as train data set in the nested 
[cell-types-eval-control-workflow](https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow-cross-validation/tree/feature/integrate_celltypes_workflow) which contains the following cell type predictors: 
* [garnett-eval-workflow](https://github.com/ebi-gene-expression-group/garnett-eval-workflow)
* [scmap-eval-workflow](https://github.com/ebi-gene-expression-group/scmap-eval-workflow) - 2 versions: 'cluster' and 'cell'
* [scpred-eval-workflow](https://github.com/ebi-gene-expression-group/scpred-eval-workflow)
* [label-analysis-eval-workflow](https://github.com/ebi-gene-expression-group/label-analysis-eval-workflow)

Workflow diagram:
![](https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow-cross-validation/blob/feature/integrate_celltypes_workflow/cross-validation-pipeline.png)

### Pipeline configuration
This workflow currently allows for the configuration of the number of folds to split the input data under the `k_folds_num` parameter.

### Triggering the pipeline
You will need Nextflow installed to run analyses. The best way to install it is via Conda. It is recommended to create a clean environment to avoid dependency conflicts. Run the following commands:

```
conda create --name nextflow
conda activate nextflow
conda install nextflow
```

In order to pull the latest version of the nested workflow containing individual pipelines you'll have to issue the following command from the `cell-types-eval-control-workflow-cross-validation` directory:
```
./bin/fetch-control-and-tool-eval-workflows.sh
```
Finally, to run the pipeline run: 
```
nextflow run main.nf --profile cluster
```
### Results
The output of the pipeline is a folds' average performance table that will be located in `data/cross_val_output/` directory. 
