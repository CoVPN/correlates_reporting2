# Cox Proportional Hazards Modeling of Correlates of Risk

## Contents 

* `code`: scripts for pre-processing and analyzing post-processed data
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results
* `output`: results files produced by statistical analyses


## Reproducibility 

The following bash scripts assume that we start at the root level of the repository.

### VaxArt mock data correlates

To generate covpn_correlates_cor_coxph_nextgen_mock_DATESTRING.pdf, run:
```{bash}
export TRIAL=nextgen_mock
cd cor_coxph
make 
```

