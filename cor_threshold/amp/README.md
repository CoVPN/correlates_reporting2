# AMP Trial NP threshold code

## Reproducibility

This project uses a project-level `renv.lock`. Setup:

- Download the repository, for example by downloading and unzipping a release from https://github.com/CoVPN/correlates_reporting2/releases.
- Use R 4.4.2.
- Place `dt703new.csv`, `dt704new.csv`, and `amp_sieve_pooled_marks_final_v9.csv` directly in this `amp/` folder.
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

The following shell command should be run from the `amp/` folder.

To render the report:

```bash
Rscript -e "rmarkdown::render('amp_threshold_report_2026Aug25.Rmd', output_file = 'amp_threshold_report_$(date +%Y%m%d).pdf')"
```

For questions, please email jpspeng@uw.edu.
