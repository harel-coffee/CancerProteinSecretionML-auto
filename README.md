# CancerProteinSecretionML
Analysis of gene expression changes in the protein secretory pathway of different cancer types using machine learning.

## Reproducing the analyses in the manuscript

### Code
The python and R code necessary to reproduce the analyses in the manuscript can be found in the [scripts](scripts) directory of this repository. View the README therein for further details on the associated scripts.

### Environments
There are two conda environment files that define the packages necessary for running the python and R scripts: `environment_python.yml` and `environment_R.yml`, respectively. Create the environments from the files using the following command:
```
conda env create -f environment_python.yml
conda env create -f environment_R.yml
```

Activate either environment using `conda activate`:
```
conda activate psp-cancer-py
```

_Note:_ The environments were built on MacOS. If you are using a different OS and experience problems when creating either of the environments, try removing the version specified after each package in the `.yml` file (e.g., change `numpy=1.18.1` to `numpy`).

### Data
Note that you will first need to retrieve the larger data files from the associated Zenodo repository prior to re-running the analyses. This is described in further detail by the README in the [data](data) directory.


## Analysis result files
The raw analysis output files can be found in the [results](results) directory.



