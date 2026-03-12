# Summary of sharing of data and code for Fong et al. “Immune correlates analysis of the ENSEMBLE single Ad26.COV2.S dose vaccine efficacy clinical trial”

**Youyi Fong, Peter Gilbert**

**November 27, 2025**

- All analyses are based on the analysis-ready data files
•	janssen_pooled_partA_data_processed_with_riskscore.csv
•	janssen_na_partA_data_processed_with_riskscore.csv
•	janssen_la_partA_data_processed_with_riskscore.csv
•	janssen_sa_partA_data_processed_with_riskscore.csv


**Code set up**

1.	Download and unzip the release, e.g.,
wget https://github.com/CoVPN/correlates_reporting2/archive/ad3decd818129b986b02b92bc33ad5576bd
2.	Modify renv.lock to R version from 4.1.2 to 4.0.4
3.	Start R at the repo root level. Install renv 0.14.0 if needed:
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.14.0.tar.gz",
  repos = NULL,
  type = "source"
)
packageVersion("renv")

4.	Run the following to install the packages via renv.
Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxxxxxxxxxxxxxx") # use your personal github access token
renv::restore() 
5.	Modify config.yml section janssen_pooled_realbAb, field data_cleaned to point the local copy of analysis-ready data file.



## Table 2, Fig 3
correlates_reporting2/cor_coxph module
Assuming we are at the repo root level, run:
```{bash}
export TRIAL=janssen_pooled_partA
cd cor_coxph
make 
```

This generates covpn_correlates_cor_coxph_janssen_pooled_partA.pdf. 
Table 3 of the manuscript includes a post-hoc analysis, which can be reproduced by this commit: https://github.com/CoVPN/correlates_reporting2/tree/96cb95dd10b86412460ff60a4dfcf1611dada59f. To make it easy to reproduce, I extracted the following R code from the commit (any package version works; replace line 3 with local data):

```{r}
library(survey)
library(kyotil)
dat_proc = read.csv('/trials/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/janssen_pooled_partA_data_processed_with_riskscore_20240226.csv')
uloqs=c(bindSpike=238.1165, pseudoneutid50=12936*0.0653) # 844.7208
for (a in c("pseudoneutid50", "bindSpike")) {
  t="Day29"
  dat_proc[[t %.% a]] <- ifelse(dat_proc[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat_proc[[t %.% a]])
}    
dat=subset(dat_proc, Trt==1 & ph1.D29==1 & Bserostatus==0)
input = twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2.D29, data=dat)
fit = svycoxph(Surv(SevereEventTimePrimaryIncludeNotMolecConfirmedD29, SevereEventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score + as.factor(Region) +  Day29pseudoneutid50 + Day29bindSpike, input)
fit2 = svycoxph(Surv(ModerateEventTimePrimaryIncludeNotMolecConfirmedD29, ModerateEventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score + as.factor(Region) +  Day29pseudoneutid50 + Day29bindSpike, input)
getFormattedSummary(list(fit, fit2), robust=T, exp=T)
```

## Fig 2
correlates_reporting2/cor_graphical module

Assume that we have run cor_coxph, which produces an Rdata file that is needed by cor_graphical.
Assuming we are at the repo root level, run 
```{bash}
export TRIAL=janssen_pooled_partA
cd cor_graphical
# all endpoints
Rscript code/cor_data_preprocess.R D29IncludeNotMolecConfirmed 
Rscript code/cor_graphics_violin_scatter.R D29IncludeNotMolecConfirmed
# severe endpoints
Rscript code/cor_data_preprocess.R D29SevereIncludeNotMolecConfirmed 
Rscript code/cor_graphics_violin_scatter.R D29SevereIncludeNotMolecConfirmed
Rscript code/cor_data_preprocess_2.R D29SevereIncludeNotMolecConfirmed
Rscript code/CoR_assay_graphics.R D29SevereIncludeNotMolecConfirmed
Rscript code/CoR_wrcdf_with_VE_lines.R D29SevereIncludeNotMolecConfirmed
# moderate endpoints
Rscript code/cor_data_preprocess.R D29ModerateIncludeNotMolecConfirmed
Rscript code/cor_graphics_violin_scatter.R D29ModerateIncludeNotMolecConfirmed
Rscript code/cor_data_preprocess_2.R D29ModerateIncludeNotMolecConfirmed
Rscript code/CoR_assay_graphics.R D29ModerateIncludeNotMolecConfirmed
Rscript code/CoR_wrcdf_with_VE_lines.R D29ModerateIncludeNotMolecConfirmed
# make pdf reports
cd ..
## overall cases
bash ./_build_chapter.sh cor_graphical 
## severe cases
### add if (COR_postfix != "IncludeNotMolecConfirmed") next to line 555 of report.rmd, comment out if (attr(config,"config")=="janssen_pooled_partA") next at line 422 & 444 & 472 and "ENSEMBLE" at line 20 of report.rmd
### edit line 61 COR_postfix_list in the report.Rmd from "IncludeNotMolecConfirmed" to "SevereIncludeNotMolecConfirmed" or "ModerateIncludeNotMolecConfirmed"
### run the following shell command and rename the generated pdf to covpn_correlates_cor_graphical_janssen_pooled_partA_severe.pdf.
bash ./_build_chapter.sh cor_graphical 
## To generate the pdf for moderate cases
### edit line 61 COR_postfix_list in the report.Rmd from "SevereIncludeNotMolecConfirmed" to "ModerateIncludeNotMolecConfirmed"
### run the following shell command and rename the generated pdf to covpn_correlates_cor_graphical_janssen_pooled_partA_moderate.pdf.
bash ./_build_chapter.sh cor_graphical 
```

Lastly, the moderate and severe-critical plots were assembled next to each other in Illustrator to make Fig 2. For D29 nAb, for instance, the figure below was assembled from Fig 2.25 in the severe report and Fig 2.25 in the moderate report.

## Table 1
correlates_reporting2/cor_tabular module

```{bash}
export TRIAL=janssen_pooled_partA
cd cor_tabular
make # Error
```
The code in this snapshot is not self-consistent. It is possible to fix it but will take some effort.


## Fig 5
correlates_reporting2/cop_exposureproximal/ensemble_severe_correlates

This project uses a project-level renv.lock. Setup:

- Download and unzip the files from https://github.com/CoVPN/correlates_reporting2/releases/tag/2.2.7

- Assume that we have R 4.4.2 installed.

- Assume that we have renv 0.13.2 installed. If not, open R console at the project level (the folder containing this readme file), and run the following commands. Note that we use renv 0.13.2, which uses renv/activate.R, instead of newer versions because of some errors with the newer versions. (If in a slurm env, load an appropriate R module and a CMmake module. The latter is needed to install some packages, e.g., nloptr, lme4.
  ```{R}
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.13.2.tar.gz",
    repos = NULL,
    type = "source"
  )
  
  packageVersion("renv")  # should show ‘0.13.2’
  ```
- Run the following R command at the project level to install package dependencies:
  ```{R}
  renv::restore()
  ```
- Look for read.csv in all R scripts and modify the code to point to the local copy of analysis-ready data file. For example, if the project is located inside the repo, use the following lines of code to read from the data file on the SCHARP file system.
  ```{R}
  config.reporting <- config::get(config = "janssen_pooled_partA", file="../../config.yml") 
  dat<-read.csv(config.reporting$data_cleaned)
  ```


To reproduce VE2_Scale_LRT2_event.pdf, run the following shell commands at the project level. The first four steps use sbatch to run jobs on a high performance cluster. 

- Run the following and wait till all the jobs finish (squeue). Estimated time minutes.
    ```{bash}
    bash ./run_step_1.sh
    ```
- Run the following and wait till all the jobs finish (squeue). Estimated time 1 hour.
    ```{bash}
    bash ./run_step_2.sh
    ```
- Run the following and wait till all the jobs finish (squeue). Estimated time 4 days, depending on the availability of the nodes.
    ```{bash}
    bash ./run_step_3.sh
    ```
- Run the following and wait till all the jobs finish. Estimated time 1 hour.
    ```{bash}
    bash ./run_step_4.sh
    ```
- Run the following:
    ```{bash}
    Rscript computeSimVE_Scale.R
    Rscript PlotFig5Revision.R
    ```

## Fig 4 , Table 4
Controlled risk and VE analysis; mediation analysis

Setup:

- Download the repository from https://github.com/Avi-Kenny/Code__IC-Pipeline/

- Assume that we have a high performance computing environment with a slurm scheduler.

- Assume that we have R 4.4.2 loaded, e.g. on the Fred Hutch cluster, ml R-bundle-CRAN/2024.11-foss-2024a

- Assume that we have additional modules required for packages installation and execution loaded, including GSL and CMake. The bundle above takes care of all modules needed.

- Assume that we have renv installed. If not, open R console at the project level (the folder containing this readme file), and run the following commands to install the current renv from CRAN, which is 1.1.5 as of September 2025. A different renv version may also work.
  ```{r}
  install.packages("renv")
  
  packageVersion("renv")  
  ```

- Run the following R command at the project level to install package dependencies. Note that superlearner and pch need to be installed separately as shown below.
  ```{R}
  renv::restore() # choose restore when presented with options
   ```


- Modify line 651 (copied below) in analysis.R to point to the local copy of analysis-ready data file.
  ```{r}
    cfg2$folder_local <- "../covpn/adata/"
  ```

Make sure that the analysis string at line 18 of analysis.R is set to "Janssen (partA)" as shown below.
```{r}
  cfg2 <- list(analysis="Janssen (partA)", calc_ests=T, seed=1)
```

To generate the controlled risk plots, run the following commands in a bash shell on the repo root level. 64 jobs are created to run in parallel on a slurm cluster. Each job renders one of the 64 Rmd files. The script run_r.sh contains the following code:
```{bash}
sbatch --array=1-64 run_r.sh
```
When all jobs are done, the controlled risk and VE plots will be in the folder Figures + Tables\Janssen (partA) plots. To get the mediation results, run the following command in R at the project level:
```{r}
source("make_mediation_table.R")
```
