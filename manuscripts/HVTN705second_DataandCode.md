# Summary of sharing of data and code for Kenny et al. “Immune correlates analysis of the Imbokodo (HVTN 705/HPX2008) efficacy trial of a mosaic HIV-1 vaccine regimen evaluated in Southern African people assigned female sex at birth: a two-phase case-control study”

**Youyi Fong, Peter Gilbert**

**July 10, 2026**

- All analyses are based on the analysis-ready data file HVTN705_secondcasecontrolprocesseddata_v12.csv unless otherwise specified.
- A copy of statistical reports are on the SCHARP network drive /trials/vaccine/p705/analysis/lab/cc/copcor/correlates_reports

## Magnitude breadth plots

/trials/vaccine/p705/analysis/lab/cc/descriptive/lum_bama_breadth_iga/code/v705_bama_breadth_iga_cc_report.Rmd

Note: code for making magnitude-breadth score and mdw score is in /trials/vaccine/p705/analysis/lab/cc/pdata/lum_bama/HVTN705_LUM05_20221223_casecontrol/code


## Fig 1 Violin plots by case status

/trials/vaccine/p705/analysis/lab/cc/descriptive/lum_bama_iga/code/v705_cc_lum_iga_pt_report.Rmd

## Fig 3A 3C

Threshold analysis.

**Statistical reports**: cor_threshold_hvtn705second.pdf The generated report looks slightly different from the figure in the paper. To exactly reproduce the paper figure, find the version of report used for the manuscript and replace the commit in the instructions below with the commit in the appendix of that report.

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/youyifong/correlates_reporting2-2/archive/4c44b7bcf9bc92a1f58d73de54df4b1177e608ad.zip
unzip 4c44b7bcf9bc92a1f58d73de54df4b1177e608ad.zip
cd correlates_reporting2-2-4c44b7bcf9bc92a1f58d73de54df4b1177e608ad

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder azd1222 points to a local copy of azd1222_stage2_data_processed_20240515.csv

# d) generate report pdf
export TRIAL=hvtn705second
cd cor_threshold
# edit Makefile to uncomment line 50
make

```


