#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#-----------------------------------------------

library(here)
library(dplyr)
library(abind)
source(here("code", "params.R"))


#load(here("..", "data_clean",paste0(attr(config, "config"), "_params.Rdata"))) # file removed. objects moved to _common.R

dat.long.twophase.sample <- readRDS(here(
  "data_clean",
  "long_twophase_data.rds"
))


MaxbAbDay29 <- ifelse(exists("MaxbAbDay29"), MaxbAbDay29, NA)
MaxbAbDay57 <- ifelse(exists("MaxbAbDay57"), MaxbAbDay57, NA)
MaxID50ID80Day29 <- ifelse(exists("MaxID50ID80Day29"), MaxID50ID80Day29, NA)
MaxID50ID80Day57 <- ifelse(exists("MaxID50ID80Day57"), MaxID50ID80Day57, NA)
Maxlive50Day29 <- ifelse(exists("Maxlive50Day29"), Maxlive50Day29, NA)
Maxlive50Day57 <- ifelse(exists("Maxlive50Day57"), Maxlive50Day57, NA)


MaxbAbB <- try(dat.long.twophase.sample %>%
  filter(assay %in% c("bindSpike", "bindRBD", "bindN",
                      "bindSpike_B.1.1.7", "bindSpike_B.1.351", "bindSpike_P.1", "bindRBD_B.1.1.7", "bindRBD_B.1.351", "bindRBD_P.1")) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(MaxbAbB) == "try-error") {  ## bAb assays are unavailable
  MaxbAbB <- NA
}

MaxID50ID80B <- try(dat.long.twophase.sample %>%
  filter(assay %in% c("pseudoneutid50", "pseudoneutid80", "pseudoneutid50sa", "pseudoneutid80la")) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(MaxID50ID80B) == "try-error") { ## nAb assays are unavailable
  MaxID50ID80B <- NA
}

Maxlive50B <- try(dat.long.twophase.sample %>%
  filter(assay %in% c("liveneutmn50")) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(Maxlive50B) == "try-error") { ## nAb assays are unavailable
  Maxlive50B <- NA
}

MaxbAb <- max(MaxbAbB, MaxbAbDay29, MaxbAbDay57, na.rm = TRUE)
MaxID50ID80 <- max(MaxID50ID80B, MaxID50ID80Day29, MaxID50ID80Day57, na.rm = TRUE)
Maxlive50 <- max(Maxlive50B, Maxlive50Day29, Maxlive50Day57, na.rm = TRUE)
# axis limits for plotting assay readouts
assay_lim <- array(NA, dim = c(length(assay_immuno), length(times), 2))
dimnames(assay_lim) <- list(assay_immuno, times, c("lb", "ub"))


assay_lim[, !grepl("Delta", times), "lb"] <- 
  floor(log10(llods[assay_immuno] / 2)) # lower bound same for all assays - days
if (study_name=="AZD1222" & grepl("bind", assays)) {
  assay_lim["bindSpike", !grepl("Delta", times), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
  }# prevent19 and AZ has llod for bAb as NA, use lloq instead
if (study_name=="PREVENT19") {
  assay_lim["bindSpike", !grepl("Delta", times), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
  assay_lim["bindRBD", !grepl("Delta", times), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
}
assay_lim[assay_immuno %in% bAb_assays, !grepl("Delta", times), "ub"] <- 
  max(ceiling(MaxbAb) + ceiling(MaxbAb) %% 2, ceiling(log10(uloqs[assay_immuno])))
assay_lim[assay_immuno %in% nAb_assays, !grepl("Delta", times), "ub"] <-
  ceiling(MaxID50ID80) + ceiling(MaxID50ID80) %% 2
assay_lim[assay_immuno %in% live_assays, !grepl("Delta", times), "ub"] <-
  ceiling(Maxlive50) + ceiling(Maxlive50) %% 2


assay_lim[, grepl("Delta", times), "lb"] <- -2 # lower bound same for all assays - delta
assay_lim[assay_immuno %in% bAb_assays, grepl("Delta", times), "ub"] <- 
  ceiling(MaxbAb - min(log10(llods[bAb_assays] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(llods[bAb_assays] / 2), na.rm=T)) %% 2

if (study_name=="AZD1222" & grepl("bind", assays)) {
  assay_lim[assay_immuno %in% bAb_assays, grepl("Delta", times), "ub"] <- 
    ceiling(MaxbAb - min(log10(lloqs[bAb_assays] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lloqs[bAb_assays] / 2), na.rm=T)) %% 2
}# prevent19 and AZ has llod for bAb as NA, use lloq instead
if (study_name=="PREVENT19") {
  assay_lim["bindSpike", grepl("Delta", times), "ub"] <- 
    ceiling(MaxbAb - min(log10(lloqs["bindSpike"] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lloqs["bindSpike"] / 2), na.rm=T)) %% 2
  
  assay_lim["bindRBD", grepl("Delta", times), "ub"] <- 
    ceiling(MaxbAb - min(log10(lloqs["bindRBD"] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lloqs["bindRBD"] / 2), na.rm=T)) %% 2
}

assay_lim[assay_immuno %in% nAb_assays, grepl("Delta", times), "ub"] <- 
  ceiling(MaxID50ID80 - min(log10(llods[nAb_assays] / 2), na.rm=T)) + ceiling(MaxID50ID80 - min(log10(llods[nAb_assays] / 2), na.rm=T)) %% 2
assay_lim[assay_immuno %in% live_assays, grepl("Delta", times), "ub"] <- 
  ceiling(Maxlive50 - min(log10(llods[live_assays] / 2))) + ceiling(Maxlive50 - min(log10(llods[live_assays] / 2))) %% 2



# Quick workaround for janssen presentation report
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  assay_lim[, !grepl("Delta", times),'lb'] <- 0
  assay_lim[, !grepl("Delta", times),'ub'] <- 3
  
  assay_lim[, grepl("Delta", times),'lb'] <- -1
  assay_lim[, grepl("Delta", times),'ub'] <- 2
}

# dup assay_lim to avoid dimention got dropped for the dataset with only one marker
if (length(assay_immuno)==1){ # i.e., AZ study
  assay_lim=abind(assay_lim, assay_lim, along=1)
}

saveRDS(assay_lim,
        file = here("data_clean", "assay_lim.rds")
)
