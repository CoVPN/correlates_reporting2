# Poisson Regression Modeling of Correlates of Risk


## Reproducibility 

This project uses the repo-level renv.lock. 


### ILiAD phase 2 correlates

Setup:
- Download and unzip the release
https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.6.zip
- Run the following R command at the repo root level to install package dependencies:
```{R}
renv::restore()
```
- Open config.yml in an editor, look for iliad_ib202p, and modify the line below to point to the local copy of analysis-ready data file.



To generate covpn_correlates_cor_poisson_iliad_ib202p.pdf, run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_poisson
make 
```

To render the R markdown files, run the following commands in a bash shell:
```{bash}
cd cor_poisson
Rscript -e "rmarkdown::render('code/comparative_immunogenicity.Rmd', output_file='../comparative_immunogenicity_$(date +%Y%m%d).pdf')"
Rscript -e "rmarkdown::render('code/diproperm_run.Rmd',              output_file='../diproperm_run_$(date +%Y%m%d).pdf')"
Rscript -e "rmarkdown::render('code/posthoc_analyses.Rmd',           output_file='../posthoc_analyses_$(date +%Y%m%d).pdf')"
```

To run Yutong's code to generate the heatmap pdf files, run the following commands in a bash shell:
```{bash}
cd cor_poisson
Rscript code/Plotting_Main_YutongJin.R
```
