#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
if (!is.null(config$assay_metadata)) {pos.cutoffs = assay_metadata$pos.cutoff}
#-----------------------------------------------

library(here)
library(dplyr)
library(abind)
source(here("code", "params.R"))

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
  filter(assay %in% bAb_assays) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(MaxbAbB) == "try-error") {  ## bAb assays are unavailable
  MaxbAbB <- NA
}

MaxID50ID80B <- try(dat.long.twophase.sample %>%
  filter(assay %in% nAb_assays) %>%
  select(B) %>%
  max(na.rm = TRUE), silent = TRUE)

if (class(MaxID50ID80B) == "try-error") { ## nAb assays are unavailable
  MaxID50ID80B <- NA
}

Maxlive50B <- try(dat.long.twophase.sample %>%
  filter(assay %in% live_assays) %>%
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
  floor(log10(lods[assay_immuno] / 2)) # lower bound same for all assays - days
if (study_name=="AZD1222" && grepl("bind", assays)) {
  assay_lim["bindSpike", !grepl("Delta", times), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
  }# prevent19 and AZ has llod for bAb as NA, use lloq instead
if (study_name=="PREVENT19") {
  assay_lim["bindSpike", !grepl("Delta", times), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
  assay_lim["bindRBD", !grepl("Delta", times), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
}
if (study_name=="VAT08") {
  assay_lim["bindSpike_mdw", !grepl("Delta", times), "lb"] <- floor(min(dat.mock[, grepl("bindSpike_mdw", colnames(dat.mock)) & !grepl("Delta", colnames(dat.mock))], na.rm=T))
  assay_lim["pseudoneutid50_mdw", !grepl("Delta", times), "lb"] <- floor(min(dat.mock[, grepl("pseudoneutid50_mdw", colnames(dat.mock)) & !grepl("Delta", colnames(dat.mock))], na.rm=T))
  for (aa in assays[grepl("bind", assays) & !grepl("mdw", assays)]){
    assay_lim[aa, !grepl("Delta", times), "lb"] <- floor(log10(lloqs[aa] / 2))
  }
}
assay_lim[assay_immuno %in% bAb_assays, !grepl("Delta", times), "ub"] <- 
  max(ceiling(MaxbAb) + ceiling(MaxbAb) %% 2, ceiling(log10(uloqs[assay_immuno][!is.infinite(uloqs[assay_immuno])])))
assay_lim[assay_immuno %in% nAb_assays, !grepl("Delta", times), "ub"] <-
  ceiling(MaxID50ID80) + ceiling(MaxID50ID80) %% 2
if (study_name=="VAT08") {
  assay_lim["pseudoneutid50_mdw", !grepl("Delta", times), "ub"] <- ceiling(max(dat.mock[,"Day22pseudoneutid50_mdw"], na.rm=T))
}
assay_lim[assay_immuno %in% live_assays, !grepl("Delta", times), "ub"] <-
  ceiling(Maxlive50) + ceiling(Maxlive50) %% 2


if (study_name=="VAT08") {assay_lim[, grepl("Delta", times), "lb"] <- -7
} else {assay_lim[, grepl("Delta", times), "lb"] <- -2 } # lower bound same for all assays - delta

assay_lim[assay_immuno %in% bAb_assays, grepl("Delta", times), "ub"] <- 
  ceiling(MaxbAb - min(log10(lods[bAb_assays] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lods[bAb_assays] / 2), na.rm=T)) %% 2

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
  ceiling(MaxID50ID80 - min(log10(lods[nAb_assays] / 2), na.rm=T)) + ceiling(MaxID50ID80 - min(log10(lods[nAb_assays] / 2), na.rm=T)) %% 2
assay_lim[assay_immuno %in% live_assays, grepl("Delta", times), "ub"] <- 
  ceiling(Maxlive50 - min(log10(lods[live_assays] / 2))) + ceiling(Maxlive50 - min(log10(lods[live_assays] / 2))) %% 2
if (study_name=="VAT08") {
  assay_lim[, grepl("Delta", times), "ub"] <- 8 
}


# Quick workaround for janssen presentation report
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  assay_lim[, !grepl("Delta", times),'lb'] <- 0
  assay_lim[, !grepl("Delta", times),'ub'] <- ifelse(attr(config,"config")=="janssen_partA_VL", 4, 3)
  
  assay_lim[, grepl("Delta", times),'lb'] <- -1
  assay_lim[, grepl("Delta", times),'ub'] <- 2
}

if(study_name=="VAT08") {# hard code for VAT08
  assay_lim["bindSpike_mdw", "B","lb"] <- -2
  assay_lim["pseudoneutid50_mdw", "B","lb"] <- -1
  assay_lim["bindSpike_mdw", "Day22","lb"] <- -5
  assay_lim["pseudoneutid50_mdw", "Day22","lb"] <- -3
  assay_lim["bindSpike_mdw", "Day43","lb"] <- -6
  assay_lim["pseudoneutid50_mdw", "Day43","lb"] <- -4
  
  assay_lim[, "Delta22overB","lb"] <- rep(-3, 15)
  assay_lim["bindSpike_mdw", "Delta22overB","lb"] <- -5
  assay_lim["pseudoneutid50_mdw", "Delta22overB","lb"] <- -4
  
  assay_lim[, "Delta43overB","lb"] <- rep(-3, 15)
  assay_lim["bindSpike_mdw", "Delta43overB","lb"] <- -5
  assay_lim["pseudoneutid50_mdw", "Delta43overB","lb"] <- -5
  
  
  assay_lim[, "B","ub"] <- 5
  assay_lim["bindSpike_mdw", "B","ub"] <- 11
  
  assay_lim[, "Day22","ub"] <- 6
  assay_lim[, "Day43","ub"] <- 5
  
  assay_lim[, "Delta22overB","ub"] <- 6
  assay_lim[, "Delta43overB","ub"] <- 4
} 

# dup assay_lim to avoid dimention got dropped for the dataset with only one marker
if (length(assay_immuno)==1){ # i.e., AZ study
  assay_lim=abind(assay_lim, assay_lim, along=1)
}

saveRDS(assay_lim,
        file = here("data_clean", "assay_lim.rds")
)
