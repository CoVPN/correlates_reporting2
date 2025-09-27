# ENSEMBLE severe COVID exposure-proximal correlate analysis


## Reproduciblity

This project uses a project-level renv.lock. Setup:

- Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases

- Assume that we have R 4.4.2 installed.

- Assume that we have renv 0.13.2 installed. If not, open R console at the project level (the folder containing this readme file), and run the following commands. Note that we use renv 0.13.2, which uses renv/activate.R, instead of newer versions because of some errors with the newer versions. (If in a slurm env, load an appropriate R module and a CMmake module. The latter is needed to install some packages, e.g., nloptr, lme4.
  ```{r}
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.13.2.tar.gz",
    repos = NULL,
    type = "source"
  )
  
  packageVersion("renv")  # should show ‘0.13.2’
  ```
- Run the following R command at the project level to install package dependencies:
  ```{R}
  renv::restore()
  ```
- Look for read.csv in all R scripts and modify the code to point to the local copy of analysis-ready data file. For example, if the project is located inside the repo, use the following lines of code to read from the data file on the SCHARP file system.
  ```{R}
  config.reporting <- config::get(config = "janssen_pooled_partA", file="../../config.yml") 
  dat<-read.csv(config.reporting$data_cleaned)
  ```


To reproduce VE2_Scale_LRT2_event.pdf, run the following shell commands at the project level. The first four steps use sbatch to run jobs on a high performance cluster. 

- Run the following and wait till all the jobs finish (squeue). Estimated time minutes.
    ```{bash}
    bash ./run_step_1.sh
    ```
- Run the following and wait till all the jobs finish (squeue). Estimated time 1 hour.
    ```{bash}
    bash ./run_step_2.sh
    ```
- Run the following and wait till all the jobs finish (squeue). Estimated time 4 days, depending on the availability of the nodes.
    ```{bash}
    bash ./run_step_3.sh
    ```
- Run the following and wait till all the jobs finish. Estimated time 1 hour.
    ```{bash}
    bash ./run_step_4.sh
    ```
- Run the following:
    ```{bash}
    Rscript computeSimVE_Scale.R
    Rscript PlotFig5Revision.R
    ```
