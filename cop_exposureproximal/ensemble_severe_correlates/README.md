# ENSEMBLE severe COVID exposure-proximal correlate analysis


## Reproduciblity

This project uses a project-level renv.lock. 

Setup:
- Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases
- Install package dependencies by opening an R console at the project level (the folder containing this readme file), and run: 
```{R}
renv::restore()
```
- Open the repo config.yml in an editor, look for iliad_ib202p, and modify the line below to point to the local copy of analysis-ready data file.


The following shell commands are to be run at the project level.

- Run the following and wait till all the jobs finish (squeue). Estimated time minutes.
    ```{bash}
    bash ./run_step_1.sh
    ```
- Run the following and wait till all the jobs finish (squeue). Estimated time 1 hour.
    ```{bash}
    bash ./run_step_2.sh
    ```
- Run the following and wait till all the jobs finish (squeue). Estimated time 1 day, depending on the availability of the nodes.
    ```{bash}
    bash ./run_step_3.sh
    ```
- Run the following:
    ```{bash}
    Rscript ComputeSimVE_Scale.R
    Rscript PlotFig5Revision.R
    ```