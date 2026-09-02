# Sanofi Stage 2 CoP manuscript NP threshold code

## Reproducibility

This project uses a project-level `renv.lock`. Setup:

- Download the repository, for example by wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.8.zip
- Use R 4.4.2.
- cd cor_threshold/sanofi_stage2
- Place `vat08_combined_data_processed_20250417.csv` directly in this folder
- Open an R console in this folder. Ignore error.
- Run the following commands to install and verify the package dependencies:

```r
renv::restore() # choose 1 
renv::status()
```

If errors occur while rendering the Rmd file, these commands can help. 
```
# install rmarkdown from prebuilt binaries. Attempts to install from source (without options line) often result in failure.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))
install.packages("rmarkdown")

install.packages("remotes")
remotes::install_github("tlverse/sl3")
remotes::install_github("tlverse/tmle3")
remotes::install_github("CoVPN/npthreshold")
install.packages("DT")
```

For maintainers, use the following commands to record package changes in `renv.lock`:

```r
renv::snapshot()
renv::status()
```

The following shell command should be run from the `sanofi_stage2/` folder.

To render the report:

```bash
Rscript -e "rmarkdown::render('sanofi_stage2_npthreshold_report_2026Aug25.Rmd', output_file = 'sanofi_stage2_npthreshold_report_$(date +%Y%m%d).pdf')"
```

For questions, please email jpspeng@uw.edu.
