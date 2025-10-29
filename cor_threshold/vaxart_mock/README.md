# VaxArt SAP mock analysis NP threshold code 

## Reproducibility

This project uses a project-level renv.lock. 

**Setup steps:**

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

4. Run the following R command at the project level to install package dependencies:
  ```{R}
  renv::restore()
  ```

**Running the code:**

1. Look for read.csv in all R/Rmd scripts and modify the code to point to the local copy of analysis-ready data file if needed.

2. To render the report html, run the following shell command at the project level.
```{bash}
Rscript -e "rmarkdown::render('cor_threshold_barda_mock.Rmd', output_file='cor_threshold_barda_mock_$(date +%Y%m%d).pdf')"
```
