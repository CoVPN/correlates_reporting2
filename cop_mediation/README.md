# Vaccine Efficacy Mediation Analyses


## Reproducibility

This project uses a module-level renv.lock. 

**Restore renv**

1. Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases

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

4. Start R in cop_mediation and run
  ```{R}
  renv::restore()
  ```

### Sanofi stage 2

1. Look for read.csv in code/mediation_report_sanofi_stage2.R and modify the code to point to the local copy of analysis-ready data file if needed.

2. Run the following shell command in cop_mediation.

```{bash}
export TRIAL=vat08_combined
export stage=2
make
```
sanofi_Day43_20250617.csv and sanofi_FR_20250617.csv are used for Table 2 in the manuscript.