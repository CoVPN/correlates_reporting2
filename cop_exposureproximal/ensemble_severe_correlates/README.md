# ENSEMBLE severe COVID exposure-proximal correlate analysis


## Reproduciblity

This project uses a project-level renv.lock. Setup:
- Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases
- Assume that we have R 4.4.2 installed.
- Assume that we have renv 0.13.2 installed. If not, open R console at the project level (the folder containing this readme file), and run the following commands. Note that we use renv 0.13.2, which uses renv/activate.R, instead of newer versions because of some errors with the newer versions. (If in a slurm env, load an appropriate R module and a CMmake module. The latter is needed to install some packages, e.g., nloptr, lme4.
)
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

