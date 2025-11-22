# Poisson Regression Modeling of Correlates of Risk



## Reproduce ILiAD phase 2 correlates results

**Setup**
- Download and unzip the release
https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.6.zip
- Run the following R command at the repo root level to install package dependencies as this project uses the repo-level renv.lock:
```{R}
renv::restore()
```
- Open config.yml in an editor, look for iliad_ib202p, and modify the line below to point to the local copy of analysis-ready data file.

**Analysis**
- This produces two pdfs: correlates_cor_poisson_iliad_ib202p_{datestring}.pdf and posthoc_analyses_{datestring}.pdf:
```{bash}
export TRIAL=iliad_ib202p
cd cor_poisson
make 
```

**Mapping between mansucript TLFs and reports**

| Manuscript | Report                        |
|------------|-------------------------------|
| Table 1    | correlates_cor_poisson Figure |
| Table 2    | posthoc_analyses Table       |
| Table 3    | Value 3                       |
| Figure 1   | Value 4                       |
| Figure 2   | Value 5                       |
| Figure 3   | Value 6                       |
| Figure 4   | Value 7                       |
