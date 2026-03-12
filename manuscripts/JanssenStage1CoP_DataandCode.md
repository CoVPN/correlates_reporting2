# Summary of sharing of data and code for Fong et al. “Immune correlates analysis of the ENSEMBLE single Ad26.COV2.S dose vaccine efficacy clinical trial”

**Youyi Fong, Peter Gilbert**

**November 27, 2025**

- All analyses are based on the analysis-ready data files: janssen_pooled_partA_data_processed_with_riskscore.csv, janssen_na_partA_data_processed_with_riskscore.csv, janssen_la_partA_data_processed_with_riskscore.csv, janssen_sa_partA_data_processed_with_riskscore.csv


**Code set up**

1.	Download and unzip the release, e.g.,
wget https://github.com/CoVPN/correlates_reporting2/archive/ad3decd818129b986b02b92bc33ad5576bd
2.	Modify renv.lock to R version from 4.1.2 to 4.0.4
3.	Start R at the repo root level. Install renv 0.14.0 if needed:
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.14.0.tar.gz",
  repos = NULL,
  type = "source"
)
packageVersion("renv")

4.	Run the following to install the packages via renv.
Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
renv::restore() 
5.	Modify config.yml section janssen_pooled_realbAb, field data_cleaned to point the local copy of analysis-ready data file.



## Figure 1, Figure 2, Table 1, and Table 2

Assuming we are at the repo root level, run:
```{bas}
export TRIAL=janssen_pooled_real
make cor_report
```

This generates _report_cor/covpn_correlates_cor_janssen_pooled_real.pdf. 
Manus Fig 1a: Report page 430 
Manus Fig 1b: Report page 432 
Manus Fig 1c: Report page 434
Manus Table 1: Report page 51 
Manus Fig 2a: Report page 818 
Manus Fig 2b: Report page 825
Manus Fig 2c: Report page 832 
Manus Table 2: Report page 812 and 813

## Figure 4

N/A. These tables and figures collate results from multiple manuscripts.


