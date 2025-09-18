# Reproducibility 

## ILiAD phase 2 correlates

To generate covpn_correlates_cor_poisson_iliad_ib202p.pdf, run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_poisson
make 
```

To render the R markdown files, run the following commands in a bash shell:
```{bash Rmd}
cd cor_poisson
Rscript -e "rmarkdown::render('code/comparative_immunogenicity.Rmd')"
Rscript -e "rmarkdown::render('code/diproperm_run.Rmd')"
Rscript -e "rmarkdown::render('code/posthoc_analyses.Rmd')"
```
