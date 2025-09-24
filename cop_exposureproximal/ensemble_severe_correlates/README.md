# ENSEMBLE severe COVID exposure-proximal correlate analysis


## Reproduciblity

This project uses a project-level renv.lock. 

Setup:
- Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases
- Install package dependencies by opening an R console at the project level (the folder containing this readme file), and run 
```{R}
renv::restore()
```
- Open the repo config.yml in an editor, look for iliad_ib202p, and modify the line below to point to the local copy of analysis-ready data file.


The following shell commands are to be run at the project level.

First, run the following command and wait till all the jobs finish (squeue).
```{bash}
bash ./run_step_1.sh
```

Second, run the following command and wait till all the jobs finish (squeue).
```{bash}
bash ./run_step_2.sh
```

