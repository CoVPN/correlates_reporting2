# Cox Proportional Hazards Modeling of Correlates of Risk

## Contents 

* `code`: scripts for pre-processing and analyzing post-processed data
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results
* `output`: results files produced by statistical analyses


## Reproducibility 

This module uses a repo-level renv.lock. Setup:
- Assume that we have R 4.0.4 installed
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
- Run the following R command at the repo level to install package dependencies:
  ```{R}
  renv::restore()
  ```


The following bash scripts assume that we start at the root level of the repository.

### VaxArt mock data correlates

To generate covpn_correlates_cor_coxph_nextgen_mock_DATESTRING.pdf, run:
```{bash}
export TRIAL=nextgen_mock
cd cor_coxph
make 
```

