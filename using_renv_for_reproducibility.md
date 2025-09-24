# Using renv for Reproducibility

renv can be used at one of three levels: repo-level, module-level, and project level. For example, cor_threshold is a module-level directory and under that, sanofi_stage2 is a project-level directory. 

## Setting up renv at the project-level

This is the most reproduciblity way of using renv because each manuscript has its own renv.lock.

To start, make sure there are three files in the project folder:

- .here  This is an empty file that is used by the here package to identify the project root.
- .gitignore  This file tells git which files/folders to ignore. It should contain the following lines:
  ```
  .html
  .pdf
  .Rhistory
  .RData
  .Rproj.user
  renv/
  ```
- readme_reproducibility.md  This file contains instructions on how to reproduce the reports.


Open a new R console in the project folder in a terminal and run the following commands. Note that we use renv 0.13.2 instead of newer versions because of familiarity with older versions (e.g., the use of renv/activate.R) and installation errors with newer versions.
```{r}
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/renv/renv_0.13.2.tar.gz",
  repos = NULL,
  type = "source"
)

packageVersion("renv")  # should show ‘0.13.2’

renv::init()
```

This will create three files/folders:

- renv/  This folder contains the local library of packages used in the project. In fact, to reduce redundancy, it contains links to the packages installed in your global R library.
- renv.lock  This is a manifest file and contains the list of packages and their versions used in the project.
- .Rprofile  This file is sourced when the project is opened and activates the renv environment It is helpful to edit the file and add a line to show the current directory so that its content is as follows:
  ```{r}
  source("renv/activate.R")
  print(getwd())
  ```
  
At this point, you can close the console.

Now open a RStudio session, change the working directory to the project folder, and click on the menu Session/Restart R. You should see a message from renv, indicating that the renv environment is activated. Alternative, if you start a new R console in the project folder, you should see the same message.

For an example of making reproducible reports, check out the project folder cor_threshold/sanofi_stage2 (https://github.com/CoVPN/correlates_reporting2/blob/master/cor_threshold/sanofi_stage2/README.md).

## Using renv at the module-level

This is not currently used.

## Using renv at the repo-level

Several modules share the same renv.lock at the repo level, including cor_graphical, cor_tabular, and cor_coxph.
