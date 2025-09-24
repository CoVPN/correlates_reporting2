# Using renv for Reproducibility


## First time using renv for a project

In the project folder, it is good to have the following three files:

- .gitignore  This file tells git which files/folders to ignore. It should contain the following lines:
```
.html
.pdf
.Rhistory
.RData
.Rproj.user
renv/
```
- README.md  This file should have a section titled Reproducibility, which details how to reproduce the reports.


If in a slurm env, load an appropriate R module and a CMmake module. The latter is needed for some packages, e.g., nloptr, lme4, WeMix

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

- renv.lock  This is a manifest file and contains the list of packages and their versions used in the project. If the project is not empty, init() looks for package dependencies, installs them, and writes the dependencies to renv.lock. If there are many package dependencies, this step can take a while.
- renv/  This folder contains the local library of packages used in the project. In fact, to reduce redundancy, it contains links to the packages installed in your global R library.
- .Rprofile  This file is sourced when the project is opened and activates the renv environment It is helpful to edit the file and add a line to show the current directory so that its content is as follows:
  ```{r}
  source("renv/activate.R")
  print(getwd())
  ```

At this point, you can close the console.
```{r}
q()
```

Now open a RStudio session, change the working directory to the project folder, and click on the menu Session/Restart R. You should see a message from renv, indicating that the renv environment is activated. Alternative, if you start a new R console in the project folder, you should see the same message.

For an example of making reproducible reports, check out the project folder cor_threshold/sanofi_stage2/README.md



## Tips for using renv

Use renv::install and not remotes or DevTools to install or update pacakges from CRAN, github, or a url, e.g.,
- renv::install("github_id/package_name")
- renv::install("survival)
- renv::install("https://cran.r-project.org/src/contrib/Archive/kyotil/kyotil_2024.7-31.tar.gz")

renv::snapshot() updates renv.lock.

If you encounter warnings like: curl: (22) The requested URL returned error: 403 rate limit exceeded. The solution is to first authenticate with your personal github access token by running the following in R:
```{r}
Sys.setenv(GITHUB_PAT = "xxxxxxxxxxxxxx")
```

If something inexplicable goes wrong, check the $HOME directory to see if there are .Rprofile and renv/. If yes, delete them and try again. a

