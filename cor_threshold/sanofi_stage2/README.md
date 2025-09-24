# Sanofi Stage 2 CoP manuscript NP threshold code 

## Reproducibility

This project uses a project-level renv.lock. 

Setup:
- Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases
- Assume that we have R 4.4.2.
- Assume that we have renv 0.13.2 installed. If not, run the following to install it:
  ```{r}
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.13.2.tar.gz",
    repos = NULL,
    type = "source"
  )
  
  packageVersion("renv")  # should show ‘0.13.2’
  ```
- Run the following R command at the repo root level to install package dependencies:
  ```{R}
  renv::restore()
  ```
  Note for admin only: to put content in renv.lock, we install package dependencies using the following commands:
  ```{r}
  renv::install("isotone")
  renv::install("tlverse/tmle3")
  renv::install("jpspeng/npthreshold")
  renv::snapshot()
  ```
- Modify the line below in sanofi_stage2_npthreshold_report_2025June25.Rmd to point to the local copy of analysis-ready data file.
  ```{r}
  df <- read.csv("vat08_combined_data_processed_20250321.csv")
  ```
  
To render the report html, run the following commands in a bash shell:
```{bash}
Rscript -e "rmarkdown::render('sanofi_stage2_npthreshold_report_2025June25.Rmd')"
```


For any questions, please email jpspeng@uw.edu. 
