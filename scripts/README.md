# Scripts

This directory contains the scripts to perform the analyses described in the manuscript and reproduce the figures.


## Re-run the ML analysis

The `runMLanalysis.py` script is used to run the machine learning analysis and generate gene scores and model scoring metric values (e.g., ROC AUC). The [python conda environment file](../environment_python.yml) necessary to run this script can be found in the top-level directory of this repository. Ensure that the environment is activated prior to running any of the scripts (`conda activate psp-cancer-py`).

*NOTE!* The script will check if results already exist in the `results/` directory of the repository, and will skip the analysis if the same result file exists. Therefore, to regenerate the results, either delete the existing output present in the `results/` directory, or set the `overwrite_results` parameter within the `runMLanalysis.py` script to `True`.


### Specify run parameters

First specify the run parameters by editing the `runMLanalysis.py` script. In the top portion of the script, the following run parameters can be specified:


#### Class variable: `ClassVar` 

Defines the variable by which samples will be grouped. Some example options for this parameter are:
- `'mutTP53'` mutations status of TP53 (p53)
- `'CancerStatus'` sample tissue type (normal vs. tumor)
- `'TumorStageMerged'` comparison between two tumor stages
- `'AllStageCombos'` same as `TumorStageMerged`, but automatically iterates through all possible pairs of tumor stages.


#### Class variable levels: `VarLevelsToKeep`

Defines the two values of the class variable to be compared.
- `['FALSE', 'TRUE']` when `ClassVar` is `'mutTP53'`
- `['Solid Tissue Normal', 'Primary solid Tumor']` when `ClassVar` is `'CancerStatus'`
- `['stage i','stage ii']` (or any pair of stages) when `ClassVar` is `'TumorStageMerged'`
- `['stage i', 'stage ii', 'stage iii', 'stage iv']` when `ClassVar` is `'TumorStageMerged'` to perform a regression analysis (instead of binary classification)

Note that when `ClassVar` is set to `AllStageCombos`, all possible tumor stage pairs will be analyzed automatically so the `VarLevelsToKeep` variable does not need to be assigned.


#### Log transform pseudocount: `logTransOffset`

Defines the pseudocount added to the gene expression (TPM) prior to log-transformation. Default value is 1.


#### Removal of low-count genes: `med_tpm_threshold`

There are a few options for handling low-count genes:
- `'none'` = don't remove any genes.
- `'zero'` = remove genes with all zeros.
- `'X%'` = where `X` is a number from 0 to 100, removes genes with median TPM in the bottom X-percentile.
- `X` = where `X` is a number, removes genes with median TPM below `X`


#### ML model scoring metric: `score_metric`
Specify the scoring metric that the model will try to maximize. See the full list of available metrics on the [scikit-learn metrics and scoring page](https://scikit-learn.org/stable/modules/model_evaluation.html).
- Examples for binary classification: `'accuracy'`, `'average_precision'`, `'f1'`, `'roc_auc'`
- Examples for regression: `'explained_variance'`, `'neg_mean_squared_error'`, `'r2'`


### Specify paths and other run settings

Some additional settings can optionally be changed in the next section of the `runMLanalysis.py` script:

#### HDF5 data file

Name of the `.h5` file containing the transcript abundance (gene TPM values) for the PSP genes across the different samples, as well as some additional sample metadata.

#### Output directory

Specify the name of the directory to which the results of the ML analysis will be written. This should not be the full directory path, just the name of the folder which exists within the top-level repository directory. The folder will be created if it does not yet exist.

#### Overwrite option

Specify whether the results should be overwritten. If `False`, the script will skip analyses for which results files already exist.

#### List of available cancer types

The TCGA cancer type abbreviations through which to iterate. Cancer types that do not have the necessary sample types or numbers to perform an analysis will automatically be skipped, so no pre-filtering needs to be done here.


### Run the analysis

The script can be run from the command prompt:
```
python runMLanalysis.py
```

## Re-run the differential expression analyses

Differential expression (DE) analyses were run in R using the `CancerDiffExpAnalysis.R` script. 

Activate the conda R environment
```
conda activate psp-cancer-r
```

Launch RStudio from the command line to ensure that the R packages from the conda environment are available for use
```
r studio &
```

From within RStudio, source the `CancerDiffExpAnalysis.R` function by either checking `Source on Save` and saving the function, or by directly entering the command
```
# Replace with the actual path to the function!
source('~/CancerProteinSecretionML/scripts/CancerDiffExpAnalysis.R')
```

Run the function for the comparison(s) of interest; for example:
```
CancerDiffExpAnalysis(cancerType='ACC', classVar='mutTP53', classVarLevels=c('FALSE', 'TRUE'), main_dir='~/CancerProteinSecretionML')
```

See the function header for details on input parameters. It is also shown below for convenience.
```
# cancerType      Character of one cancer type, or list of multiple cancer types to analyze.
#                 Specifying 'all' will analyze all available cancer types in the data.
#                 If NULL, the function will return a list of all available cancer types.
#
# classVar        Class variable by which to separate the samples. For example:
#                 'CancerStatus', 'TumorStage', 'TumorStageMerged', or 'mutTP53'
#                 If NULL, the function will return a list of all available class variables.
#
# classVarLevels  The TWO values of the class variable by which the samples will be grouped.
#                 Note that the fold-changes will be calculated as L2/L1.
#                 For example, if classVar is 'CancerStatus', classVarLevels would likely be
#                 c('Solid Tissue Normal', 'Primary solid Tumor').
#                 If NULL, the function will return a list of all available levels for the 
#                 chosen class variable.
#
# gene            Leave empty (NULL) to perform a differential expression analysis.
#                 If a gene name is specified, then the expression of that gene in the
#                 two groups will be visualized with a plot.
#
# main_dir        Path to the CancerProteinSecretionML directory.
#                 e.g., 'usr/yourname/Documents/CancerProteinSecretionML'
```

## Re-generate the manuscript figures

The figures were generated in R using the `generate_figures.R` script. 

Activate the conda R environment
```
conda activate psp-cancer-r
```

Launch RStudio from the command line to ensure that the R packages from the conda environment are available for use
```
r studio &
```

Edit the top of the script to specify relevant paths on your local machine
```
# specify directory information
proj_dir <- '~/CancerProteinSecretionML'
fig_dir <- file.path(proj_dir, 'doc', 'manuscript', 'figures', 'fig_pieces')
results_folder <- 'results'  # do not specify full path, it is assumed to be in proj_dir
```

Either run the entire script to generate all figures, or run pieces of the script for individual figures (after first running the `Load and organize data` section).






