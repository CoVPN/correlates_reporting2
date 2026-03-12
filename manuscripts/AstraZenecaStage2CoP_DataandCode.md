# Summary of sharing of data and code for Janes et al. “Correlates of Severe and Delta COVID-19 in a Phase 3 Trial of the AZD1222 Vaccine”

**Youyi Fong, Peter Gilbert**

**March 10, 2026**

- All analyses are based on the analysis-ready data file azd1222_stage2_data_processed_20240515.csv unless otherwise specified.
- A copy of statistical reports are on the SCHARP network drive T:\covpn\p3002\analysis\correlates\stage2\reports\



## Table 1
Summary statistics of binding and neutralizing antibody markers at D57 among severe COVID-19 cases, Delta COVID-19 cases, and controls in the Day 57 analysis set.

**Statistical reports**: covpn_correlates_cor_tabular_azd1222_stage2.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/last_release_in_R40.zip
unzip last_release_in_R40.zip
cd correlates_reporting2-last_release_in_R40

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder azd1222 points to a local copy of azd1222_stage2_data_processed_20240515.csv

# d) generate report pdf
export TRIAL=azd1222_stage2
cd cor_tabular
make

```



## Figure 1 and Figure 2

**Statistical reports**: covpn_correlates_cor_graphical_azd1222_stage2.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/yiwenlu/correlates_reporting2/archive/3b85affab2cf25c78b8f6b507420d8237625129e.zip
unzip 3b85affab2cf25c78b8f6b507420d8237625129e.zip
cd correlates_reporting2-3b85affab2cf25c78b8f6b507420d8237625129e

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()
    
    renv::install('wCorr')

# c) edit config.yml so that the data_cleaned field uder azd1222 points to a local copy of azd1222_stage2_data_processed_20240515.csv

# d) generate report pdf
export TRIAL=azd1222_stage2
cd cor_graphical
make

```



## Figure 3 and Figure 4

**Statistical reports**: covpn_correlates_cor_coxph_azd1222_stage2.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/last_release_in_R40.zip
unzip last_release_in_R40.zip
cd correlates_reporting2-last_release_in_R40

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz") # if needed
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder azd1222 points to a local copy of azd1222_stage2_data_processed_20240515.csv

# d) generate report pdf
export TRIAL=azd1222_stage2
cd cor_coxph

# modify code/cor_coxph_azd1222_stage2.R
  # change line 137 to:
  # dat.plac = NULL,

make 

```