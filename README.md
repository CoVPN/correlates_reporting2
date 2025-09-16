# Generalized Correlates Analysis Reporting

## Summary

This repository houses modular workflows for statistical analyses of correlates
of risk / protection and automated reporting of analytic results. It serves as
a generalized suite of tools, based on the analyses originally designed for the
USG Biostatistics Response Team's analysis of COVID-19 vaccine efficacy trials
(archived
[here](https://github.com/CoVPN/correlates_reporting_usgcove_archive/)). See
below for brief descriptions of each of the analysis modules. This repository is
designed as the second part of an analytic pipeline, with the [correlates
processing](https://github.com/CoVPN/correlates_processing) module serving as an
upstream component.

[![Build Status](https://app.travis-ci.com/CoVPN/correlates_reporting2.svg?branch=master)](https://app.travis-ci.com/CoVPN/correlates_reporting2)
_Note:_ automated builds of the correlates of risk and protection analyses are
evaluated by the [Travis CI](https://travis-ci.org/) continuous integration
service and the PDF reports posted to this repository's [`gh-pages`
branch](https://github.com/CoVPN/correlates_reporting2/tree/gh-pages).

## List of Analysis Modules

* Correlates of Risk (CoR) Analyses
  * `cor_coxph`: Cox proportional hazards modeling of risk.
  * `cor_tabular`: Tabular descriptions of correlates of risk.
  * `cor_graphical`: Graphical descriptions of correlates of risk.
  * `cor_threshold`: Risk modeling based on correlate thresholds.
  * `cor_nonlinear`: Nonlinear modeling and evaluation.
  * `cor_surrogates`: Optimal surrogates analyses.
* Correlates of Protection (CoP) Analyses
  * `cop_prinstrat`: Principal stratification analyses.
  * `cop_stochastic`: Stochastic risk and vaccine efficacy evaluation.
  * `cop_mediation`: Correlate-mediated vaccine efficacy and risk.


## Instructions for Use

* All analysis code are written in R and we use renv to manage package versions.
* After cloning the repo, start R in the root directory. Enter 
```r
renv::restore()
```
to install package dependencies. Installation may take a few hours depending on internet speed. If renv errors occur, check to make sure that under the home directory there is no .Rprofile or libs.
* To run the analyses of a specific module on a dataset, enter the module, specify the TRIAL label corresponding to the dataset, and run Make, e.g.
```bash
cd cor_coxph
export TRIAL=janssen_pooled_partA
make
```
The dataset corresponding to janssen_pooled_partA can be found in the config.yml file at the repo root. 



## Reproducibility Guide for Investigators Contributing Modules

* Portability: One should be able to move the code to a new location and run it. Consider using [`here`](https://here.r-lib.org/) to help achieve protability.

* The location of the analysis data file should be in a single place in the whole module, most typically in a config.yml file (used through the R config package).

* The versions of R and packages used should be managed with [`renv`](https://rstudio.github.io/renv/).

* Either a Makefile or a bash script should be included so that the analyses and report generation can be run with a single command. 


---

## License

The contents of this repository are distributed under the GPL-3 license. See
file [`LICENSE.md`](https://github.com/CoVPN/correlates_reporting2/blob/master/LICENSE.md)
for details.
