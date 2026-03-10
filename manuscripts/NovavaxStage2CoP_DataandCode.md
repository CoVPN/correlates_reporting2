# Summary of sharing of data and code for Fong et al. “Analysis of Antibody Markers as Immune Correlates of Risk of Severe COVID-19 in the PREVENT-19 Efficacy Trial of the NVX-CoV2373 Recombinant Protein Vaccine”

**Youyi Fong, Peter Gilbert**

**March 10, 2026**

- All analyses are based on the analysis-ready data file prevent19_stage2_data_processed_20250325.csv unless otherwise specified.
- A copy of statistical reports are on the SCHARP network drive T:\covpn\p3005\analysis\correlates\stage2\reports\



## Table 1
List of Symptoms Defining the Severe COVID-19 Endpoint.

N/A



## Table 2
D35 Delta and Ancestral Antibody Marker Positive Response Frequencies and Geometric Means by COVID-19 Outcome Status.

**Statistical reports**: covpn_correlates_cor_tabular_prevent19_stage2.pdf

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

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_stage2_data_processed_20250325.csv

# d) generate report pdf
export TRIAL=prevent19_stage2
cd cor_tabular
make

```


## Figure 1
Violin plots showing D35 levels.

**Statistical reports**: covpn_correlates_cor_graphical_prevent19_stage2.pdf

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

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_stage2_data_processed_20250325.csv

# d) generate report pdf
export TRIAL=prevent19_stage2
cd cor_graphical

# modify report.Rmd
  # change line 40-41 to:
  # times_ = c("Day35")
  # labels.time = c("Day 35")

make

```



## Figure 2, Figure 3, and Table 3

**Statistical reports**: covpn_correlates_cor_coxph_prevent19_stage2.pdf

To **reproduce** this report, follow the following instructions:

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

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_stage2_data_processed_20250325.csv

# d) generate report pdf
export TRIAL=prevent19_stage2
cd cor_coxph

# modify code/cor_coxph_prevent19_stage2.R
  # change line 160 to:
  # dat.plac = NULL,
# modify report.Rmd
  # change line 10 to:
  # show.q=F; has.scaled=T; has.alt=F; show.sample.size=T; show.forestplots=F; show.ve.curves=F; special.note=""; has.plac=F; plot.geq=F; show.tertile.curves=T
  # change line 20 to:
  # assays = c("pseudoneutid50_D614G", "pseudoneutid50_Delta", "bindSpike_D614", "bindSpike_Delta1") 


make 

```

## Figure 4
Hazard ratio with reference antibody level LOD/2 (for nAb-ID50) or LLOQ/2 (for Spike IgG) of Delta COVID-19 as a function of current antibody marker level modeled over time.

**Statistical reports**: PreventProximalCorrelateNew.pdf

To reproduce this report, follow the following instructions on https://github.com/yinghuang124/Exposure-Proximal
