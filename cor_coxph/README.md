# Cox Proportional Hazards Modeling of Correlates of Risk

## Contents 

* `code`: scripts for pre-processing and analyzing post-processed data
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results
* `output`: results files produced by statistical analyses


## Reproducibility 

This module uses a repo-level renv.lock. 

**Setup steps:**

1. Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases

2. Assume that we have R 4.0.4 installed (or loaded on a high performance cluster).

3. Assume that we have renv 0.13.2 installed. If not, open R console at the project level (the folder containing this readme file) and run the following command.
)
  ```{r}
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.13.2.tar.gz",
    repos = NULL,
    type = "source"
  )
  
  packageVersion("renv")  # should show ‘0.13.2’
  ```

4. Run the following R command at the project level to install package dependencies:
  ```{R}
  renv::restore()
  ```


**Running the code:**

1. Look for read.csv in all R/Rmd scripts and modify the code to point to the local copy of analysis-ready data file if needed.

2. To generate the report pdf, run the following command from the root level of the repository.

**VaxArt mock data correlates**

```{bash}
export TRIAL=nextgen_mock
cd cor_coxph
make 
```

**ENSEMBLE severe correlates**
```{bash}
export TRIAL=janssen_pooled_partA
cd cor_coxph
make 
