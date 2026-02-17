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


### ENSEMBLE severe correlates
```{bash}
export TRIAL=janssen_pooled_partA
cd cor_coxph
make 
```


### Sanofi Stage 2 correlates manuscript
```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/sanofi_stage2_correlates_R40.zip
unzip sanofi_stage2_correlates_R40.zip
cd correlates_reporting2-sanofi_stage2_correlates_R40

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field under vat08_combined points to a local copy of
#    vat08_combined_data_processed_20250417.csv

# d) generate report pdf
export TRIAL=vat08_combined
export stage=2    
cd cor_coxph
make 
```


### COV2008 T cell correlates manuscript
```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/last_release_in_R40.zip
unzip last_release_in_R40.zip
cd correlates_reporting2-last_release_in_R40

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field under cov2008_tcell points to a local copy of 
#    cov2008_tcell_data_processed_20251211.csv

# d) generate report pdf
export TRIAL=cov2008_tcell
cd cor_coxph
make 
```


### COVAIL T cell correlates manuscript
```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/last_release_in_R40.zip
unzip last_release_in_R40.zip
cd correlates_reporting2-last_release_in_R40

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field under covail_tcell points to a local copy of
#    covail_data_processed_20250818.csv

# d) generate report pdf
export TRIAL=covail_tcell
cd cor_coxph
make 
```
