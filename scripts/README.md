# Scripts

This directory contains the scripts to perform the analyses described in the manuscript and reproduce the figures.


## Re-run the ML analysis

The `runMLanalysis.py` script is used to re-run the machine learning analysis and generate gene scores and ROC AUC values. The [conda environment file](../environment.yml) necessary to run this script can be found in the top-level directory of this repository.

*NOTE!* The script will check if results already exist in the `results/` directory of the repository, and will skip the analysis if the same result file exists. Therefore, to regenerate the results, first delete the existing output present in the `results/` directory.


### Specify run parameters

First specify the run parameters by editing the `runMLanalysis.py` script. In the top portion of the script, specify the following parameters:


#### Class variable: `ClassVar` 

Defines the variable by which samples will be grouped. Some example options for this parameter are:
- `'mutTP53'` mutations status of TP53 (p53)
- `'CancerStatus'` sample tissue type (normal vs. tumor)
- `'TumorStageMerged'` comparison between two tumor stages
- `'AllStageCombos'` same as `TumorStageMerged`, but automatically iterates through all possible pairs of tumor stages.


#### Class variable levels: `VarLevelsToKeep`

Defines the two values of the class variable to be compared.
- `['FALSE', 'TRUE']` when `ClassVar` is `'mutTP53'`  or `'Mutations'`
- `['Solid Tissue Normal', 'Primary solid Tumor']` when `ClassVar` is `'CancerStatus'`
- `['stage i','stage ii']` (or any pair of stages) when `ClassVar` is `'TumorStageMerged'`

Note that when `ClassVar` is set to `AllStageCombos`, all possible tumor stage pairs will be analyzed automatically so the `VarLevelsToKeep` variable does not need to be assigned.


#### Log transform pseudocount: `logTransOffset`

Defines the pseudocount added to the gene expression (TPM) prior to log-transformation. Default value is 1.


#### Removal of low-count genes: `med_tpm_threshold`

There are a few options for handling low-count genes:
- `'none'` = don't remove any genes.
- `'zero'` = remove genes with all zeros.
- `'X%'` = where `X` is a number from 0 to 100, removes genes with median TPM in the bottom X-percentile.
- `X` = where `X` is a number, removes genes with median TPM below `X`


### Run the analysis

The script can be run from the command prompt:
```
python runMLanalysis.py
```

## Re-run the differential expression analyses


## Re-generate the manuscript figures






