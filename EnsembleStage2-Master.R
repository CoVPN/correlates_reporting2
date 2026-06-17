require(renv)
# renv::load()
dir <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/code"
renv::restore(dir)

require(lubridate)
require(survminer)
require(CFsurvival)
require(knitr)
require(sas7bdat)
require(survival)
require(MASS)
require(rms)
require(haven)
require(epitools)
require(ranger)
require(vioplot)
require(latex2exp)
require(mvtnorm)
require(hal9001)
require(CFsurvival)

set.seed(206)
spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
dir.create(paste0(spath, "/reports"), showWarnings = FALSE, recursive=TRUE)


source(paste0(spath, "/code/Methods-CausalCumulativeIncidence.R"))
source(paste0(spath, "/code/Methods-AverageTreatmentEffect.R"))
source(paste0(spath, "/code/Methods-NPCorrelates.R"))

source(paste0(spath, "/code/Demographics.R"))
source(paste0(spath, "/code/Demographics-immuno.R"))

source(paste0(spath, "/code/Endpoints-adj.R"))
source(paste0(spath, "/code/Endpoints-adj-LatinAmerica.R"))
source(paste0(spath, "/code/Endpoints-adj-US.R"))
source(paste0(spath, "/code/Endpoints-adj-SouthAfrica.R"))
source(paste0(spath, "/code/Endpoints-Immuno-adj.R"))

source(paste0(spath, "/code/HybridEfficacy.R"))
source(paste0(spath, "/code/HybridEfficacy-US.R"))
source(paste0(spath, "/code/HybridEfficacy-LatinAmerica.R"))
source(paste0(spath, "/code/HybridEfficacy-SouthAfrica.R"))

source(paste0(spath, "/code/Immunogenicity.R"))

source(paste0(spath, "/code/Correlates-MarkerDist.R"))
source(paste0(spath, "/code/Correlates-MarkerDist-US.R"))
source(paste0(spath, "/code/Correlates-MarkerDist-LatinAmerica.R"))
source(paste0(spath, "/code/Correlates-MarkerDist-RSA.R"))

source(paste0(spath, "/code/Correlates-CoxRegression.R"))
source(paste0(spath, "/code/Correlates-CoxRegression-US.R"))
source(paste0(spath, "/code/Correlates-CoxRegression-LatinAmerica.R"))
# insufficient number of endpoints to run Cox regression in RSA
# source(paste0(spath, "/code/Correlates-CoxRegression-RSA.R"))

source(paste0(spath, "/code/Correlates-NPRegression.R"))
