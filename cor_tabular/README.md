#  Tabular Descriptions of Correlates of Risk



## Reproducibility 

This module uses the repo-level renv.lock.

## ILiAD phase 2 correlates

```{bash}
export TRIAL=iliad_ib202p
cd cor_tabular
make
```



### Sanofi stage2 correlates

```{bash}
# a) obtaining the code
wget https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/xxxxxxxxx
unzip xxxxx
cd xxxxx

# b) restore R package dependencies
R
    Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
    renv::restore()

# c) edit config.yml so that the data_cleaned field uder vat08_combined points to a local copy of vat08_combined_data_processed_20250417.csv

# d) generate report pdf
export TRIAL=vat08_combinedv
cd cor_tabular
make
```
