# Reproducibility

## ILiAD phase 2 correlates

To generate covpn_correlates_cor_tabular_iliad_ib202p.pdf, run the following commands in a bash shell:
```{bash}
export TRIAL=iliad_ib202p
cd cor_tabular
make
```

To render the R markdown files, run the following commands in a bash shell:
```{bash Rmd}
bash _build_chapter.sh cor_tabular
```

