# ENSEMBLE severe COVID exposure-proximal correlate analysis

## Reproduciblity

This project uses a project-level renv.lock. 

**Setup steps:**

1. Download and unzip the files from https://github.com/CoVPN/correlates_reporting2/releases/tag/2.2.7

2. Assume that we have R 4.4.2 installed (or loaded on a high performance cluster).

3. Assume that we have renv 1.1.5 installed. If not, open R console at the project level (the folder containing this readme file) and run the following command.
)
  ```{r}
  # if the following does not work, it is because 1.1.5 is the current release. Then try install.packages("renv")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz",
    repos = NULL,
    type = "source"
  )
  
  packageVersion("renv")  
  ```

4. Run the following R command at the project level to install package dependencies:
  ```{R}
  renv::restore()
  ```

**Running the code:**

1. Look for read.csv in all R/Rmd scripts and modify the code to point to the local copy of analysis-ready data file if needed.

2. To reproduce VE2_Scale_LRT2_event.pdf, run the following shell commands at the project level. The first four steps use sbatch to run jobs on a high performance cluster. 

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

