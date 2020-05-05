#!/usr/bin/env bash
CROSSVAL_WORKFLOW_ROOT="$PWD"
CONTROL_EVAL_WORKFLOW="$CROSSVAL_WORKFLOW_ROOT/cell-types-eval-control-workflow"
#EVAL_WORKFLOWS="$CONTROL_EVAL_WORKFLOW/cell-types-eval-workflows"

#Nested submodule structure where eval-workflow containing submodules for individual pipelines is located within the control workflow and this last into the outer cross-validation workflow

# Clone or update the cell-types-eval-control-workflow repo 
if [ ! -d 'cell-types-eval-control-workflow' ]; then
   git clone --recursive https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow.git $CONTROL_EVAL_WORKFLOWS
fi

pushd $CONTROL_EVAL_WORKFLOWS > /dev/null
git checkout feature/cross_validation > /dev/null
git pull origin feature/cross_validation > /dev/null
git submodule update --recursive --remote > /dev/null
popd > /dev/null
