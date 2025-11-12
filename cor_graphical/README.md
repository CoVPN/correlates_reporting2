#  Graphical Descriptions of Correlates of Risk

## Contents

* `code`: scripts for pre-processing and analyzing post-processed data
* `data_clean`: post-processed data used as input to statistical analyses
* `figs`: visualizations of exploratory and analytic results, 
          1. violin & line plots, intercurrent vs per protocol, case vs non-case, at (Day 1), Day 29 Day 57
          2. scatter plots, marker vs age, intercurrent vs per protocol, case vs non-case, at (Day 1) Day 29 Day 57
* `output`: results files produced by statistical analyses
* `slurm`: scheduler scripts for submission of batch jobs


## Reproducibility 

This module uses the repo-level renv.lock.

### ILiAD phase 2 correlates

There are four reports. 

- To generate covpn_correlates_cor_graphical_iliad_ib202p.pdf, replace the "C0iliad_ib202p_PPAI_C14only" in the lines 27 and 28 of the ./Makefile and the line 91 of the ./report.Rmd to be "C0iliad_ib202p", run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_graphical
make
```

- To generate covpn_correlates_cor_graphical_iliad_ib202p.pdf, replace the "C0iliad_ib202p_PPAI_C14only" in the lines 27 and 28 of the ./Makefile and the line 91 of the ./report.Rmd to be "C0iliad_ib202p_C14only", run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_graphical
make
```

- To generate covpn_correlates_cor_graphical_iliad_ib202p.pdf, replace the "C0iliad_ib202p_PPAI_C14only" in the lines 27 and 28 of the ./Makefile and the line 91 of the ./report.Rmd to be "C0iliad_ib202p_PPAI", run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_graphical
make
```

- To generate covpn_correlates_cor_graphical_iliad_ib202p.pdf, make sure the it is "C0iliad_ib202p_PPAI_C14only" that is included in the lines 27 and 28 of the ./Makefile and the line 91 of the ./report.Rmd, run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_graphical
make
```

### VaxArt correlates

- To generate covpn_correlates_cor_graphical_nextgen_mock.pdf, run the following commands in a bash shell:

for short report, make sure line 62 of the report.Rmd shows "nextgen_output_flag = 1"
```{bash}
export TRIAL=nextgen_mock
cd cor_graphical
make
```

for full report, make sure line 62 of the report.Rmd shows "nextgen_output_flag = 3"
```{bash}
export TRIAL=nextgen_mock
cd cor_graphical
make
```
