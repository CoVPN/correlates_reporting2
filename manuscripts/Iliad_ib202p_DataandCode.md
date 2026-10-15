# Reproduce ILiAD phase 2 correlates results

**Setup**
- Use R 4.4

- Download and unzip the release
https://github.com/CoVPN/correlates_reporting2/archive/refs/tags/2.2.9.zip
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

Manuscript	Report
Table 1	    posthoc_analyses_20261015.pdf, Table 2.1
Table 2	    posthoc_analyses_20261015.pdf, Table 5.2
Table 3	    posthoc_analyses_20261015.pdf, Table under Figure 6.1 (p. 31)
Figure 1	  comparative_immunogenicity_20260514.pdf, Figure 4.1-4.4
Figure 2	  posthoc_analyses_20261015.pdf, Figure 5.3
Figure 3	  correlates_cor_poisson_iliad_ib202p_20260512.pdf: A) Figure 1.2; B) Figure 1.6. 
Figure 4	  posthoc_analyses_20261015.pdf, Figure 6.1
