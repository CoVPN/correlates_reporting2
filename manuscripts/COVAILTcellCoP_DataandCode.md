# Summary of sharing of data and code for Zhang et al. “Distinct T-cell correlates of COVID-19 risk before and after a second booster in previously SARS-CoV-2 infected vs. naïve COVAIL trial participants ”

**Bo Zhang, Youyi Fong, Bhavesh Borate, Lindsay Carpp, Peter Gilbert**

**March 12, 2026**

- All analyses are based on the analysis-ready data file covail_data_processed_20250818.csv unless otherwise specified.
- A copy of statistical reports are on the SCHARP network drive T:\covpn\COVAILcorrelates\analysis\correlates\reports\McElrathTcellsarticle




## Figure 1, Figure 2
Distributions of T-cell responses at D1 and at D15, shown separately for each one-dose vaccine second booster arm, in a representative random cohort (non-cases plus a subset of cases) of the naïve cohort. 

**Statistical reports**: covpn_correlates_cor_voilin_covail_tcell_09042025.pdf

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

# c) edit config.yml so that the data_cleaned field uder covail_tcell points to a local copy of covail_tcell_data_processed_20250818.csv

# d) generate report pdf
export TRIAL=covail_tcell
cd cor_graphical
make
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_graphical_$TRIAL.pdf', config_file = '_bookdown_cor_graphical.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Table 1

CD4+ and CD8+ T-cell responses among the one-dose mRNA group, for COVID-19 endpoint cases (booster-proximal and booster-distal pooled) and non-cases and for the naïve and non-naïve cohorts. 

**Statistical reports**: covpn_correlates_cor_tabular_covail_tcell.pdf

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

# c) edit config.yml so that the data_cleaned field uder covail_tcell points to a local copy of covail_tcell_data_processed_20250818.csv

# d) generate report pdf
export TRIAL=covail_tcell
cd cor_tabular
Rscript code/make_table_all.R D35
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_tabular_$TRIAL.pdf', config_file = '_bookdown_cor_tabular.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Table 2
Superlearner results on best models predicting overall COVID-19 (booster-proximal and booster-distal cases pooled) in the one-dose mRNA group, based on neutralizing antibody and T-cell responses, for (A) the naïve cohort and (B) the non-naïve cohort.

**Statistical reports**: covail_tcell_multivariable_SuperLearner_reports

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

# c) edit config.yml so that the data_cleaned field uder covail_tcell points to a local copy of covail_tcell_data_processed_20250818.csv

# d) generate report pdf
export TRIAL=covail_tcell
cd cor_surrogates
make

```




## Figure 2 

Scatterplots of pairs ofD35 antibodymarker values (spike IgG, RBDIgG, pseudovirus nAb-ID50) for baseline SARS-CoV-2 negative per-protocol vaccine recipients in the immunogenicity subcohort (U.S. study sites).

**Statistical reports**: covpn_correlates_immuno_graphical_covail_tcell.pdf

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
    
# c) edit config.yml so that the data_cleaned field uder covail_tcell points to a local copy of covail_tcell_data_processed_20250818.csv

# d) generate report pdf
R
    # two additional packages are needed for immuno_graphical
    # dummies needs to installed separately from an archive
    renv::install('https://cran.r-project.org/src/contrib/Archive/dummies/dummies_1.5.6.tar.gz')
    renv::install('SWIM)

# edit code/descriptive_graphics_two_phase_plots.R to comment out line 23, which fails to run

export TRIAL=covail_tcell
cd immuno_graphical
make
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_immuno_graphical_$TRIAL.pdf', config_file = '_bookdown_immuno_graphical.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```


## Figure 3
Covariate-adjusted hazard ratios (HRs) 

**Statistical reports**: covpn_correlates_cor_coxph_covail_tcell_20250729.pdf

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

# c) edit config.yml so that the data_cleaned field uder covail_tcell points to a local copy of covail_tcell_data_processed_20250818.csv

# d) generate report pdf
export TRIAL=covail_tcell
cd cor_coxph
make 

```

## Figure 4 
Controlled risk (covariate-adjusted cumulative incidence) of COVID-19 by T-cell response level at D1 or at D15 in the one-dose mRNA group, for the naïve and non-naïve cohorts

**Statistical reports**: covpn_correlates_controlledrisk_covail_06132025_McElrathTcells.pdf

To reproduce this report, follow the following instructions. 

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

# c) edit config.yml so that the data_cleaned field uder covail_tcell points to a local copy of covail_tcell_data_processed_20250818.csv

# d) generate report pdf
export TRIAL=covail_tcell
cd cop_contrisk

```



## Figure 5
T-cell subset expression patterns at D1 and at D15 in booster-proximal cases vs. non-cases in the one-dose mRNA group, for the naïve and non-naïve cohorts

**Statistical reports**: COMPASS/cyt_combination.pdf, LeidenClustering



## Figure 6
Bivariate analysis of nAb-ID50 BA.4/5 and CD4+ T-cell IFN- and/or IL-2 Spike BA.4/5 responses at D1 or D15 as correlates of risk in the one-dose mRNA group

**Statistical reports**: 