# Summary of sharing of data and code for Luedtke et al. “Immune Correlates Analysis of Antibody Responses Against SARS-CoV-2 Variants in the ENSEMBLE Vaccine Efficacy Trial”

**Youyi Fong, Peter Gilbert**

**November 27, 2025**

- All analyses are based on the analysis-ready data file janssen_partA_VL_data_processed_20240226.csv


**Code set up**


1.	Download and unzip the release
https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.3.zip
2.	Start R at the repo root level. Run the following to install the packages via renv.
```{R}
Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
renv::restore() 
renv::install('youyifong/kyotil') # update to current version
renv::install('CoVPN/copcor') # update to current version
renv::install('fmsb')
renv::install('wCorr')
renv::install('ggnewscale')
```

3.	Modify config.yml line 254 to point the local copy of analysis-ready data file.




## Figure 2, 3 and 4

**Statistical report**: covpn_correlates_cor_coxph_janssen_partA_VL.pdf

To reproduce the report, in cor_coxph/report_trial/report_ensemble.Rmd, line 118, add
```{R}
      show.tertile.curves <- show.ve.curves <- has.plac <- T; plot.geq <- F
```
In cor_coxph/code/cor_coxph_janssen_partA_VL.R, change line 48 to:
```{r}
    write(paste0(gsub("_", " ", a, fixed = TRUE),     " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", a, ".txt"))
```

Run:
```{bash}
cd cor_coxph
export TRIAL=janssen_partA_VL
make
```


## Figure 1

**Statistical report**: covpn_correlates_cor_graphical_janssen_partA_VL.pdf

To reproduce the report, run:
```{bash}
cd cor_graphical
export TRIAL=janssen_partA_VL
make
cd ..
bash _build_chapter.sh cor_graphical
```


## Table 1

**Statistical report**: covpn_correlates_cor_tabular_janssen_partA_VL.pdf

To reproduce the report, code from correlates_reporting2-2.2.7 is needed
```{bash}
cd cor_graphical
export TRIAL=janssen_partA_VL
make
cd ..
bash _build_chapter.sh cor_tabular
```

