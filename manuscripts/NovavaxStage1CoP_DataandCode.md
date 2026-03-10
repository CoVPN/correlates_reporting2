# Summary of sharing of data and code for Fong et al. “Immune correlates analysis of the PREVENT-19 COVID-19 vaccine efficacy clinical trial”

**Youyi Fong, Peter Gilbert**

**March 6, 2026**

- All analyses use code in this release if not otherwise specified: https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.tar.gz
- All analyses are based on the analysis-ready data file prevent19_data_processed_20221016.csv unless otherwise specified.
- A copy of statistical reports are on the SCHARP network drive T:\covpn\p3005\analysis\correlates\Part_A_Blinded_Phase_Data\reports\




## Figure 1
D35 antibody marker level by COVID-19 outcome status in baseline SARS-CoV-2 negative per-protocol vaccine recipients (U.S. study sites).

**Statistical reports**: covpn_correlates_cor_graphical_prevent19.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_data_processed_20221016.csv

# d) generate report pdf
export TRIAL=prevent19
cd cor_graphical
make
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_graphical_$TRIAL.pdf', config_file = '_bookdown_cor_graphical.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Table 1

D35 antibody marker SARS-CoV-2 seroresponse rates and geometric means in the U.S. cohort by COVID-19 outcome status

**Statistical reports**: covpn_correlates_cor_tabular_prevent19.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_data_processed_20221016.csv

# d) generate report pdf
export TRIAL=prevent19
cd cor_tabular
Rscript code/make_table_all.R D35
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_tabular_$TRIAL.pdf', config_file = '_bookdown_cor_tabular.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Figure 2 

Scatterplots of pairs ofD35 antibodymarker values (spike IgG, RBDIgG, pseudovirus nAb-ID50) for baseline SARS-CoV-2 negative per-protocol vaccine recipients in the immunogenicity subcohort (U.S. study sites).

**Statistical reports**: covpn_correlates_immuno_graphical_prevent19.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()
    
# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_data_processed_20221016.csv

# d) generate report pdf
R
    # two additional packages are needed for immuno_graphical
    # dummies needs to installed separately from an archive
    renv::install('https://cran.r-project.org/src/contrib/Archive/dummies/dummies_1.5.6.tar.gz')
    renv::install('SWIM)

# edit code/descriptive_graphics_two_phase_plots.R to comment out line 23, which fails to run

export TRIAL=prevent19
cd immuno_graphical
make
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_immuno_graphical_$TRIAL.pdf', config_file = '_bookdown_immuno_graphical.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```


## Figure 3, Table 2, and Figure 5

**Statistical reports**: covpn_correlates_cor_coxph_prevent19.pdf

To **reproduce** this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_data_processed_20221016.csv

# d) generate report pdf
export TRIAL=prevent19
cd cor_coxph
make 
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_coxph_$TRIAL.pdf', config_file = '_bookdown_cor_coxph.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```

## Figure 4 

Analyses of D35 antibodymarkers as a correlate of risk in baseline SARSCoV-2 negative per-protocol vaccine recipients (U.S. study sites)

**Statistical reports**: covpn_correlates_cor_threshold_prevent19.pdf

To reproduce this report, follow the following instructions. Notes: 1) There are extra instructions in step d). 2) The generated report pdf is located in _report_cor. 

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder prevent19 points to a local copy of prevent19_data_processed_20221016.csv

# d) generate report pdf
export TRIAL=prevent19
cd cor_threshold
Rscript code/params.R D35
Rscript code/clean_data.R D35
Rscript code/Run_Threshold_analysis.R D35

# modify _common.R 
  # change line 995 to: labels <- c(" ", labels)
# change cor_threshold/report.md:
  # replace TRIAL <- "janssen_trial_real" with TRIAL <- Sys.getenv("TRIAL")

Rscript code/plotting.R D35
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_threshold_$TRIAL.pdf', config_file = '_bookdown_cor_threshold.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Table 3 and Figure 6

N/A. These tables and figures collate results from multiple manuscripts.
