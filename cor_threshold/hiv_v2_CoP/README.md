# HIV V2 Antibody CoP manuscript NP threshold code

This folder contains code for the non-parametric threshold analysis in the HIV V2 Antibody correlates of protection manuscript.

## Reproducibility

This project uses a project-level `renv.lock`. Setup:

- Download the repository, for example by downloading and unzipping a release from https://github.com/CoVPN/correlates_reporting2/releases.
- Use R 4.4.2.
- Place `controlledVEdata702.csv`, `HVTN705_secondcasecontrolprocesseddata_v14.csv`, `controlledVEdata505_lab.csv`, and `controlledVEdataRV144_v2.csv` directly in this `hiv_v2_CoP/` folder.
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

The following shell command should be run from the `hiv_v2_CoP/` folder.

To render the report:

```bash
Rscript -e "rmarkdown::render('threshold_analysis_hiv_2026Aug25.Rmd', output_file = 'threshold_analysis_hiv_$(date +%Y%m%d).html')"
```

For questions, please email jpspeng@uw.edu.
