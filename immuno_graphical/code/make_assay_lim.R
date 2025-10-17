#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
if (!is.null(config$assay_metadata)) {pos.cutoffs = assay_metadata$pos.cutoff}
if (study_name == "VaxArt_Mock") {
  assays = assays[!grepl("nasal|IgG_nasal|saliva|bindN_IgA", assays)]
  assay_metadata = assay_metadata %>% filter(!grepl("nasal|saliva|bindN_IgA", assay))
} # will remove later
#-----------------------------------------------

library(here)
library(dplyr)
library(abind)
source(here::here("code", "params.R"))

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
assay_lim <- array(NA, dim = c(length(assay_immuno), length(times_), 2))
dimnames(assay_lim) <- list(assay_immuno, times_, c("lb", "ub"))


assay_lim[, !grepl("Delta", times_), "lb"] <- 
  floor(log10(lods[assay_immuno] / 2)) # lower bound same for all assays - days
if (study_name=="AZD1222" && grepl("bind", assays)) {
  assay_lim["bindSpike", !grepl("Delta", times_), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
  }# prevent19 and AZ has llod for bAb as NA, use lloq instead
if (attr(config,"config")=="prevent19") {
  assay_lim["bindSpike", !grepl("Delta", times_), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
  assay_lim["bindRBD", !grepl("Delta", times_), "lb"] <- floor(log10(lloqs["bindSpike"] / 2))
}
if (attr(config,"config")=="prevent19_stage2") {
  assay_lim["pseudoneutid50_D614G", , "lb"] <- floor(log10(lods["pseudoneutid50_D614G"] / 2))
  assay_lim["pseudoneutid50_Delta", , "lb"] <- floor(log10(lods["pseudoneutid50_Delta"] / 2))
  assay_lim["bindSpike_D614", , "lb"] <- floor(log10(lloqs["bindSpike_D614"] / 2))
  assay_lim["bindSpike_Beta", , "lb"] <- floor(log10(lloqs["bindSpike_Beta"] / 2))
  assay_lim["bindSpike_Alpha", , "lb"] <- floor(log10(lloqs["bindSpike_Alpha"] / 2))
  assay_lim["bindSpike_Gamma", , "lb"] <- floor(log10(lloqs["bindSpike_Gamma"] / 2))  
  assay_lim["bindSpike_Delta1", , "lb"] <- floor(log10(lloqs["bindSpike_Delta1"] / 2))
  assay_lim["bindSpike_Delta3", , "lb"] <- floor(log10(lloqs["bindSpike_Delta3"] / 2))
  assay_lim["bindSpike_Delta2", , "lb"] <- floor(log10(lloqs["bindSpike_Delta2"] / 2))
  assay_lim["bindSpike_Omicron", , "lb"] <- floor(log10(lloqs["bindSpike_Omicron"] / 2))
}
if (study_name=="VAT08") {
  assay_lim["bindSpike_mdw", !grepl("Delta", times_), "lb"] <- 0
  assay_lim["pseudoneutid50_mdw", !grepl("Delta", times_), "lb"] <- 0
  for (aa in assays[grepl("bind", assays) & !grepl("mdw", assays)]){
    assay_lim[aa, !grepl("Delta", times_), "lb"] <- floor(log10(lloqs[aa] / 2))
  }
  for (aa in assays[grepl("pseudo", assays) & !grepl("mdw", assays)]){
    assay_lim[aa, !grepl("Delta", times_), "lb"] <- floor(log10(lods[aa] / 2))
  }
}
assay_lim[assay_immuno %in% bAb_assays, !grepl("Delta", times_), "ub"] <- 
  max(ceiling(MaxbAb) + ceiling(MaxbAb) %% 2, ceiling(log10(uloqs[assay_immuno][!is.infinite(uloqs[assay_immuno])])))
assay_lim[assay_immuno %in% nAb_assays, !grepl("Delta", times_), "ub"] <-
  ceiling(MaxID50ID80) + ceiling(MaxID50ID80) %% 2
#if (study_name=="VAT08") {
#  assay_lim["pseudoneutid50_mdw", !grepl("Delta", times_), "ub"] <- ceiling(max(dat_proc[,"Day22pseudoneutid50_mdw"], na.rm=T))
#}
assay_lim[assay_immuno %in% live_assays, !grepl("Delta", times_), "ub"] <-
  ceiling(Maxlive50) + ceiling(Maxlive50) %% 2


if (study_name=="VAT08") {assay_lim[, grepl("Delta", times_), "lb"] <- -3.5
} else {assay_lim[, grepl("Delta", times_), "lb"] <- -2 } # lower bound same for all assays - delta

assay_lim[assay_immuno %in% bAb_assays, grepl("Delta", times_), "ub"] <- 
  ceiling(MaxbAb - min(log10(lods[bAb_assays] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lods[bAb_assays] / 2), na.rm=T)) %% 2

if (study_name=="AZD1222" && grepl("bind", assays)) {
  assay_lim[assay_immuno %in% bAb_assays, grepl("Delta", times_), "ub"] <- 
    ceiling(MaxbAb - min(log10(lloqs[bAb_assays] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lloqs[bAb_assays] / 2), na.rm=T)) %% 2
}# prevent19 and AZ has llod for bAb as NA, use lloq instead
if (attr(config,"config")=="prevent19") {
  assay_lim["bindSpike", grepl("Delta", times_), "ub"] <- 
    ceiling(MaxbAb - min(log10(lloqs["bindSpike"] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lloqs["bindSpike"] / 2), na.rm=T)) %% 2
  
  assay_lim["bindRBD", grepl("Delta", times_), "ub"] <- 
    ceiling(MaxbAb - min(log10(lloqs["bindRBD"] / 2), na.rm=T)) + ceiling(MaxbAb - min(log10(lloqs["bindRBD"] / 2), na.rm=T)) %% 2
}

assay_lim[assay_immuno %in% nAb_assays, grepl("Delta", times_), "ub"] <- 
  ceiling(MaxID50ID80 - min(log10(lods[nAb_assays] / 2), na.rm=T)) + ceiling(MaxID50ID80 - min(log10(lods[nAb_assays] / 2), na.rm=T)) %% 2
assay_lim[assay_immuno %in% live_assays, grepl("Delta", times_), "ub"] <- 
  ceiling(Maxlive50 - min(log10(lods[live_assays] / 2))) + ceiling(Maxlive50 - min(log10(lods[live_assays] / 2))) %% 2
if (study_name=="VAT08") {
  assay_lim[, grepl("Delta", times_), "ub"] <- 4 
}


# Quick workaround for janssen presentation report
if(study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") {
  assay_lim[, !grepl("Delta", times_),'lb'] <- 0
  assay_lim[, !grepl("Delta", times_),'ub'] <- ifelse(attr(config,"config")=="janssen_partA_VL", 4, 3)
  
  assay_lim[, grepl("Delta", times_),'lb'] <- -1
  assay_lim[, grepl("Delta", times_),'ub'] <- 2
}

if(study_name=="VAT08" & F) {# hard code for VAT08
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

if (attr(config,"config")=="prevent19_stage2") {
  assay_lim[, ,"ub"] = 5
}

if (attr(config,"config")=="nextgen_mock") {
  assay_lim[, , "ub"] = 7
  assay_lim[, , "lb"] = 2
  
  assay_lim[c("pseudoneutid50_sera", "pseudoneutid50_sera_XBB.1.5", "pseudoneutid50_sera_KP.2", "pseudoneutid50_sera_mdw"), , "ub"] = 6
  
  assay_lim[c("T4_IFNg_OR_IL2_N_Index", "T4_IFNg_OR_IL2_Spike_KP.2", "T4_IFNg_OR_IL2_Spike_Index.D614", "T8_IFNg_OR_IL2_N_Index", "T8_IFNg_OR_IL2_Spike_KP.2", "T8_IFNg_OR_IL2_Spike_Index.D614"), c("B", "Day31", "Day91", "Day181", "Day366"), "ub"] = 1
  
  assay_lim[c("T4_IFNg_OR_IL2_N_Index", "T4_IFNg_OR_IL2_Spike_KP.2", "T4_IFNg_OR_IL2_Spike_Index.D614", "T8_IFNg_OR_IL2_N_Index", "T8_IFNg_OR_IL2_Spike_KP.2", "T8_IFNg_OR_IL2_Spike_Index.D614"), c("Delta31overB", "Delta91overB", "Delta181overB", "Delta366overB"), "ub"] = 1
  
  assay_lim[c("T4_IFNg_OR_IL2_N_Index", "T4_IFNg_OR_IL2_Spike_KP.2", "T4_IFNg_OR_IL2_Spike_Index.D614", "T8_IFNg_OR_IL2_N_Index", "T8_IFNg_OR_IL2_Spike_KP.2", "T8_IFNg_OR_IL2_Spike_Index.D614"), c("B", "Day31", "Delta31overB"), "lb"] = -2.6
  
  assay_lim[c("T4_IFNg_OR_IL2_N_Index", "T4_IFNg_OR_IL2_Spike_KP.2", "T4_IFNg_OR_IL2_Spike_Index.D614", "T8_IFNg_OR_IL2_N_Index", "T8_IFNg_OR_IL2_Spike_KP.2", "T8_IFNg_OR_IL2_Spike_Index.D614"), c("Delta91overB", "Delta181overB", "Delta366overB"), "lb"] = -1

  assay_lim[c("T4_IFNg_OR_IL2_N_Index", "T4_IFNg_OR_IL2_Spike_KP.2", "T4_IFNg_OR_IL2_Spike_Index.D614", "T8_IFNg_OR_IL2_N_Index", "T8_IFNg_OR_IL2_Spike_KP.2", "T8_IFNg_OR_IL2_Spike_Index.D614"), c("Day91", "Day181", "Day366"), "lb"] = 0
}

# dup assay_lim to avoid dimention got dropped for the dataset with only one marker
if (length(assay_immuno)==1){ # i.e., AZ study
  assay_lim=abind(assay_lim, assay_lim, along=1)
}

saveRDS(assay_lim,
        file = here("data_clean", "assay_lim.rds")
)
