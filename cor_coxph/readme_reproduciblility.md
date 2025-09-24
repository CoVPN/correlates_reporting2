# Reproducibility 

The following bash scripts assume that we start at the root level of the repository.

## VaxArt mock data correlates

To generate covpn_correlates_cor_coxph_nextgen_mock_DATESTRING.pdf, run:
```{bash}
export TRIAL=nextgen_mock
cd cor_coxph
make 
```

