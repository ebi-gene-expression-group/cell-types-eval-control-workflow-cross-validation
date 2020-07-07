#!/usr/bin/env bash
CROSSVAL_WORKFLOW_ROOT="$PWD"
CONTROL_EVAL_WORKFLOW="$CROSSVAL_WORKFLOW_ROOT/cell-types-eval-control-workflow"
EVAL_WORKFLOWS="$CONTROL_EVAL_WORKFLOW/cell-types-eval-workflows"

# Clone or update the cell-types-eval-control-workflow repo 
if [ ! -d 'cell-types-eval-control-workflow' ]; then
   git clone --recursive https://github.com/ebi-gene-expression-group/cell-types-eval-control-workflow.git $CONTROL_EVAL_WORKFLOWS
fi

pushd $CONTROL_EVAL_WORKFLOWS > /dev/null
git checkout develop > /dev/null
git pull origin develop > /dev/null
git submodule update --recursive --remote > /dev/null
popd > /dev/null

# Clone or update the cell-types-eval-workflows repo containing submodules for individual pipelines
if [ ! -d 'cell-types-eval-workflows' ]; then
    git clone --recursive https://github.com/ebi-gene-expression-group/cell-types-eval-workflows $EVAL_WORKFLOWS
fi

pushd $EVAL_WORKFLOWS > /dev/null
git checkout origin/develop > /dev/null
git pull origin develop > /dev/null
git submodule update --recursive --remote > /dev/null
popd > /dev/null
