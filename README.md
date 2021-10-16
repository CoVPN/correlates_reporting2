# Generalized Correlates Analysis Reporting

## Summary

This repository houses modular workflows for statistical analyses of correlates
of risk and protection, and the automated reporting of analytic results. It is
a generalized suite of tools based on the analyses originally developed for the
analysis of the COVID-19 vaccine efficacy trials (archived at
https://github.com/CoVPN/correlates_reporting_usgcove_archive/). See below for
brief descriptions of each of the analysis modules.

## Contents

* Correlates of Risk (CoR) Analyses
  * `cor_tabular`: Tabular descriptions of correlates of risk.
  * `cor_graphical`: Graphical descriptions of correlates of risk.
  * `cor_coxph`: Cox proportional hazards modeling of risk.
  * `cor_threshold`: Risk modeling based on correlate thresholds.
  * `cor_nonlinear`: Nonlinear modeling and evaluation.
  * `cor_surrogates`: Optimal surrogates analyses.
* Correlates of Protection (CoP) Analyses
  * `cop_prinstrat`: Principal stratification analyses.
  * `cop_stochastic`: Stochastic risk and vaccine efficacy evaluation.
  * `cop_mediation`: Correlate-mediated vaccine efficacy and risk.

## Collaboration Guide

* [Code style guide](https://style.tidyverse.org/), with some modifications;
  this will largely be enforcd with [`styler`](https://styler.r-lib.org/).
* Project organization: _mostly_ independent subdirectories, each incorporating
  [`here`](https://here.r-lib.org/) for path resolution.
* Package version control and virtual environments using
  [`renv`](https://rstudio.github.io/renv/).
* Code review procedure: see our [contribution
   guidelines](https://github.com/CoVPN/correlates_reporting2/blob/master/CONTRIBUTING.md).

---

## License

The contents of this repository are distributed under the GPL-3 license. See
file [`LICENSE.md`](https://github.com/CoVPN/correlates_reporting2/blob/master/LICENSE.md)
for details.
