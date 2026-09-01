# Sanofi Stage 2 CoP manuscript NP threshold code

## Reproducibility

This project uses a project-level `renv.lock`. Setup:

- Download the repository, for example by downloading and unzipping a release from https://github.com/CoVPN/correlates_reporting2/releases.
- Use R 4.4.2.
- Place `vat08_combined_data_processed_20250417.csv` directly in this `sanofi_stage2/` folder.
- Open an R console in this folder.
- Run the following commands to install and verify the package dependencies:

```r
renv::restore()
renv::status()
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
