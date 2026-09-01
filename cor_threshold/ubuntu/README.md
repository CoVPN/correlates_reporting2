# UBUNTU Trial NP threshold code

## Reproducibility

This project uses a project-level `renv.lock`. Setup:

- Download the repository, for example by downloading and unzipping a release from https://github.com/CoVPN/correlates_reporting2/releases.
- Use R 4.4.2.
- Place `COVID_realdata_ics_bama_nab_bcp20241210_processed_with_riskscore.csv`, `dt_COMPASS_Spike BA4-5.csv`, and `dt_COMPASS_Spike Ancestral.csv` directly in this `ubuntu/` folder.
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

The following shell command should be run from the `ubuntu/` folder.

To render the report:

```bash
Rscript -e "rmarkdown::render('ubuntu_threshold_analysis_2026June24.Rmd', output_file = 'ubuntu_threshold_analysis_$(date +%Y%m%d).pdf')"
```

For questions, please email jpspeng@uw.edu.
