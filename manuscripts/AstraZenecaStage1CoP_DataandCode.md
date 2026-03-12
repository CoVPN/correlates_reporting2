# Summary of sharing of data and code for Benkeser et al. “Immune correlates analysis of a phase 3 trial of the AZD1222 (ChAdOx1 nCoV-19) vaccine”

**Youyi Fong, Peter Gilbert**

**March 10, 2026**

- All analyses are based on the analysis-ready data files azd1222_data_processed_20220425.csv and azd1222_bAb_data_processed_20220506.csv unless otherwise specified.
- A copy of statistical reports are on the SCHARP network drive T:\covpn\p3002\analysis\correlates\Part_A_Blinded_Phase_Data\reports\




## Figure 1
D57 antibody marker level by COVID-19 outcome status in baseline SARS-CoV-2 negative vaccine recipients.

**Statistical reports**: covpn_correlates_cor_graphical_azd1222_bAb.pdf and covpn_correlates_cor_graphical_azd1222.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()

# c) edit config.yml so that the data_cleaned field under azd1222 points to a local copy of azd1222_data_processed_20220425.csv and the data_cleaned field under azd1222_bAb points to a local copy of azd1222_bAb_data_processed_20220506.csv.

# d) generate report pdf
export TRIAL=azd1222
cd cor_graphical
make
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_graphical_$TRIAL.pdf', config_file = '_bookdown_cor_graphical.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

export TRIAL=azd1222_bAb
cd cor_graphical
make
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_graphical_$TRIAL.pdf', config_file = '_bookdown_cor_graphical.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Table 1
D57 antibody markera SARS-CoV-2 seroresponse rates and geometric means by COVID-19 outcome status.

**Statistical reports**: covpn_correlates_cor_tabular_azd1222.pdf and covpn_correlates_cor_tabular_azd1222_bAb.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()

# c) edit config.yml so that the data_cleaned field under azd1222 points to a local copy of azd1222_data_processed_20220425.csv and the data_cleaned field under azd1222_bAb points to a local copy of azd1222_bAb_data_processed_20220506.csv.

# d) generate report pdf
export TRIAL=azd1222
cd cor_tabular
Rscript code/make_table_all.R D29
Rscript code/make_table_all.R D57
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_tabular_$TRIAL.pdf', config_file = '_bookdown_cor_tabular.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

export TRIAL=azd1222_bAb
cd cor_tabular
Rscript code/make_table_all.R D29
Rscript code/make_table_all.R D57
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_tabular_$TRIAL.pdf', config_file = '_bookdown_cor_tabular.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Figure 2 
Scatterplot of D57 spike IgG vs. D57 pseudovirus (PsV)-nAb ID50 values for baseline SARS-CoV-2 negative vaccine recipients in the immunogenicity subcohort.

**Statistical reports**: pairs_Day57_Markers_BaselineNeg_vaccine_arm_AZD1222.png

This is a one-off correlation plot. To reproduce this report, follow the following instructions:



## Figure 3, Table 2, Figure 4(a)-(b), and Figure 5

**Statistical reports**: covpn_correlates_cor_coxph_azd1222_20220513.pdf and covpn_correlates_cor_coxph_bAb_azd1222_20220513.pdf

To reproduce this report, follow the following instructions:

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()

# c) edit config.yml so that the data_cleaned field under azd1222 points to a local copy of azd1222_data_processed_20220425.csv and the data_cleaned field under azd1222_bAb points to a local copy of azd1222_bAb_data_processed_20220506.csv.

# d) generate report pdf
export TRIAL=azd1222
cd cor_coxph

# modify ../config.yml
  # remove line 322:   multivariate_assays: [bindSpike+pseudoneutid50]
# modify code/cor_coxph.R
  # change line 147 to:
  # if(F) {
# modify report_by_COR.Rmd
  # change line 69 to:
  # if (F) {
  
make 
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_coxph_$TRIAL.pdf', config_file = '_bookdown_cor_coxph.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"


export TRIAL=azd1222_bAb
cd cor_coxph
make 
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_coxph_$TRIAL.pdf', config_file = '_bookdown_cor_coxph.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```

## Figure 4(c)-(d)


**Statistical reports**: covpn_correlates_cor_npthreshold_azd1222_nAb.pdf and covpn_correlates_cor_npthreshold_azd1222_nbAb.pdf

To reproduce this report, follow the following instructions. Note that the generated report pdf is located in _report_cor. 

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.0.zip
unzip 2.2.0.zip
cd correlates_reporting2-2.2.0

# b) restore R package dependencies
# Assume that we have R 4.0.4 installed
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    # Assume that we have renv 0.14.0. installed. If not, run the next line
    # install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz")
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder azd1222 points to a local copy of azd1222_data_processed_20221208.csv

# d) generate report pdf
export TRIAL=azd1222
cd cor_threshold
Rscript code/params.R D29
Rscript code/clean_data.R D29
Rscript code/Run_Threshold_analysis.R D29
Rscript code/params.R D57
Rscript code/clean_data.R D57
Rscript code/Run_Threshold_analysis.R D57

# modify _common.R 
  # change line 995 to: labels <- c(" ", labels)
# change cor_threshold/report.md:
  # replace TRIAL <- "janssen_trial_real" with TRIAL <- Sys.getenv("TRIAL")
  # add line 92: keys= c("D29", "D57") 

Rscript code/plotting.R D29
Rscript code/plotting.R D57
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_threshold_$TRIAL.pdf', config_file = '_bookdown_cor_threshold.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"


export TRIAL=azd1222_bAb
cd cor_threshold
Rscript code/params.R D29
Rscript code/clean_data.R D29
Rscript code/Run_Threshold_analysis.R D29
Rscript code/params.R D57
Rscript code/clean_data.R D57
Rscript code/Run_Threshold_analysis.R D57

Rscript code/plotting.R D29
Rscript code/plotting.R D57
cd ..
Rscript -e "bookdown::render_book(input = 'index_cor.Rmd', output_file = 'covpn_correlates_cor_threshold_$TRIAL.pdf', config_file = '_bookdown_cor_threshold.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

```



## Table 3 and Figure 6

N/A. These tables and figures collate results from multiple manuscripts.
