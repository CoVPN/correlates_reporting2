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
