# VaxArt Mock Cumulative Incidence Curves

## Reproducibility

This project uses a project-level renv.lock. Setup:

- Download the repository, e.g. download and unzip a release from https://github.com/CoVPN/correlates_reporting2/releases

- Assume that we have R 4.4.2 loaded.

- Assume that we have renv 1.1.5 installed. If not, open R console at the project level (the folder containing this readme file) and run the following command.
)
  ```{R}
  # if the following does not work, it is because 1.1.5 is the current release. Then try install.packages("renv")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz",
    repos = NULL,
    type = "source"
  )
  
  packageVersion("renv")  
  ```

- Run the following R command at the project level to install package dependencies:
  ```{R}
  renv::restore()
  ```

To render the report, run the following commands in a bash shell:
```{bash}
Rscript -e "rmarkdown::render('cumulative_incidence_report_vaxart_mock.Rmd', output_file='cumulative_incidence_report_vaxart_mock_$(date +%Y%m%d).pdf')"
```


For any questions, please email bzhang3@fredhutch.org. 
