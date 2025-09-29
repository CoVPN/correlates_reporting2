# Sanofi Stage 2 CoP manuscript NP threshold code 

## Reproducibility

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

- Modify the line below in sanofi_stage2_npthreshold_report_2025June25.Rmd to point to the local copy of analysis-ready data file.
  ```{r}
  df <- read.csv("vat08_combined_data_processed_20250321.csv")
  ```

- Note for admin only: see https://github.com/CoVPN/correlates_reporting2/blob/master/using_renv_for_reproducibility.md for more info about setting up renv for a new project. 
  
The following shell commands are to be run at the project level.

To render the report html, run the following commands in a bash shell:
```{bash}
Rscript -e "rmarkdown::render('sanofi_stage2_npthreshold_report_$(date +%Y%m%d).Rmd')"
```


For any questions, please email jpspeng@uw.edu. 
