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
| Table 1    | posthoc_analyses_20260103.pdf, Table 5.2 |
| Table 2    | posthoc_analyses_20260103.pdf, Table 2.1 |
| Table 3    | posthoc_analyses_20260103.pdf, Table 5.9, 5.14 |
| Figure 1   | comparative_immunogenicity_20260103.pdf, Figure 4.1-4.4 |
| Figure 2   | posthoc_analyses_20260103.pdf, Figure 5.3 |
| Figure 3   | correlates_cor_poisson_iliad_ib202p_20260103.pdf, Figure 1.6, 1.2 |
| Figure 4   | posthoc_analyses_20260103.pdf, Figure 6.1 |
