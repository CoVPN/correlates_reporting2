#-----------------------------------------------
# obligatory to append to the top of each script
# renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
library(cowplot) # for function plot_grid
library(grid)
library(gridExtra)
library(wCorr) # weighted correlation
#library(ggnewscale) # for new_scale_color() 
library(survey)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
library(tidyverse)


Sys.setenv(TRIAL="vat08_combined")
COR="D22D43omi"

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("code", "process_violin_pair_functions.R"))
Sys.setenv(DESCRIPTIVE = 1)
source(here::here("..", "_common.R"))
source(here::here("code", "params.R")) # load parameters

# for the order of figure panels
if ("order_in_panel" %in% colnames(assay_metadata)){
    assay_order = assay_metadata %>% dplyr::arrange(panel, order_in_panel) %>% select(assay_label_short) %>% pull()
    assay_metadata = assay_metadata %>%
        mutate(assay_label_short = factor(assay_label_short,
                                          levels = assay_order
        ))
}

dat.longer.immuno.subset.plot1 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.rds"))
if(attr(config,"config") == "vat08_combined"){
    dat.longer.immuno.subset.plot1_stage1_stage2 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1_stage1_stage2.rds"))}
if(attr(config,"config")=="janssen_partA_VL") {
    dat.longer.immuno.subset.plot1.2 <- readRDS(here::here("data_clean", "longer_immuno_data_plot1.2.rds"))}
dat.immuno.subset.plot3 <- readRDS(here::here("data_clean", "twophase_data.rds")); dat.immuno.subset.plot3$all_one <- 1 # as a placeholder for strata values
dat.immuno.subset.plot3 <- dat.immuno.subset.plot3 %>%
    mutate(Trt = factor(Trt, levels = c(0, 1), labels = c("Placebo", "Vaccine")),
           nnaive = factor(Bserostatus,
                            levels = c(0, 1),
                            labels = bstatus.labels))


get_gm <- function(dat, v, subs, sub.by, strata, weights, subset){
    rgm <- NULL
    dat_twophase <- dat %>% 
        group_by_at(strata) %>% 
        mutate(ph1cnt=n(), ph2cnt=sum(!!as.name(subset), na.rm = T)) %>% 
        filter(ph1cnt!=0 & ph2cnt!=0) %>% 
        select_at(gsub("`", "", c("Ptid", strata, weights, subset, sub.by, v, subs)))
    
    design.full <- twophase(id=list(~Ptid, ~Ptid), 
                            strata=list(NULL, as.formula(sprintf("~%s", strata))),
                            weights=list(NULL, as.formula(sprintf("~%s", weights))),
                            method="simple",
                            subset=as.formula(sprintf("~%s", subset)),
                            data=dat_twophase)
    for (i in v){
        design.ij <- subset(design.full, eval(parse(text=sprintf("!is.na(%s)", i))))
        for (j in subs){
            # cat(i,"--",j,"\n")
            ret <- svyby(as.formula(sprintf("~%s", i)),
                         by=as.formula(sprintf("~%s", paste(c(j, sub.by), collapse="+"))),
                         design=design.ij,
                         svymean, vartype="ci", na.rm=T)
            
            retn <- dat %>%
                dplyr::filter(!!as.name(subset) & !is.na(!!as.name(i))) %>%
                group_by_at(gsub("`", "", c(j, sub.by))) %>%
                summarise(N=n(), .groups="drop")
            
            rgm <- bind_rows( 
                inner_join(ret, retn, by = gsub("`", "", c(j, sub.by))) %>% 
                    rename(mag=!!as.name(i), Group=!!as.name(j)) %>% 
                    mutate(subgroup=!!j, mag_cat=!!i,
                           `GMT/GMC`= sprintf("%.2f\n(%.2f, %.2f)", 10^mag, 10^ci_l, 10^ci_u)),
                rgm)
            
        }
    }
    # rgm <- inner_join(rgm, distinct(labels_all, mag_cat, Visit, Marker), by = "mag_cat") 
    return(rgm)
}

# ###### Set 2 plots: Longitudinal plots D1, D22, D43
for (panel in c("pseudoneutid50", "bindSpike")){
    # by naive/non-naive, vaccine/placebo
    
    if (attr(config,"config") %in% c("janssen_partA_VL", "janssen_pooled_partA")) next # janssen_partA_VL doesn't need these plots
    
    if (attr(config,"config")=="vat08_combined" & panel=="pseudoneutid50") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.nAb==1)
    } else if (attr(config,"config")=="vat08_combined" & panel=="bindSpike") {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1 %>% filter(ph2.immuno.bAb==1)
    } else {
        dat.longer.immuno.subset.plot1_ = dat.longer.immuno.subset.plot1
    }
}

set2_assays = assays
if (attr(config,"config")=="vat08_combined" & dat.longer.immuno.subset.plot1_$Trialstage[1]==1) {
    # longitudinal plots for stage 1 and stage 2, non-naive ppt

    
    # Generate a full table with all the estimates: GMTR
    # non-naive ppts only
    # Vaccine vs. Placebo, by Stage 1 and 2
    
    sub.by <- c("all")
    ds.i <- filter(dat_proc %>%
                       mutate(all = 1) %>% # fake column to satisfy sub.by
                       filter(Bserostatus == 1) %>% # only do this for non-naive ppts
                       mutate(Trt=ifelse(Trt==1, "Vaccine", "Placebo")), ph1.immuno==1)
    gm.v <- apply(expand.grid(times_[!grepl("Delta",times_)], assays), 1, paste0, collapse="")
    
    subs <- "Trt"
    comp.i <- c("Vaccine", "Placebo")
    
    
    rgm_stage1_bAb_vac <- get_gm(ds.i %>% filter(Trialstage == 1), 
                                    gm.v[grepl("bindSpike", gm.v)], subs, sub.by, strata="Wstratum", 
                                    weights="wt.immuno.bAb", subset="ph2.immuno.bAb") %>% mutate(Trialstage = "Stage 1")
    rgm_stage1_nAb_vac <- get_gm(ds.i %>% filter(Trialstage == 1), 
                                    gm.v[grepl("pseudoneutid50", gm.v)], subs, sub.by, strata="Wstratum", 
                                    weights="wt.immuno.nAb", subset="ph2.immuno.nAb") %>% mutate(Trialstage = "Stage 1") 
    rgm_stage2_bAb_vac <- get_gm(ds.i %>% filter(Trialstage == 2), 
                                    gm.v[grepl("bindSpike", gm.v)], subs, sub.by, strata="Wstratum", 
                                    weights="wt.immuno.bAb", subset="ph2.immuno.bAb") %>% mutate(Trialstage = "Stage 2") 
    rgm_stage2_nAb_vac <- get_gm(ds.i %>% filter(Trialstage == 2), 
                                    gm.v[grepl("pseudoneutid50", gm.v)], subs, sub.by, strata="Wstratum", 
                                    weights="wt.immuno.nAb", subset="ph2.immuno.nAb") %>% mutate(Trialstage = "Stage 2")
    
    nab_stg1 <- dat_proc %>% filter(Trialstage == 1, ph2.immuno.nAb, Trt==1, Bserostatus == 0) %>% mutate(wtD22=Day22pseudoneutid50_10*wt.immuno.nAb, wtD43=Day43pseudoneutid50_10*wt.immuno.nAb)
    
    bab_stg1 <- dat_proc %>% filter(Trialstage == 1, ph2.immuno.bAb, Trt==1, Bserostatus == 0) %>% mutate(wtD22=Day22bindSpike*wt.immuno.bAb, wtD43=Day43bindSpike*wt.immuno.bAb)
    
    
    for (asy in set2_assays){
        f_2 <-
            f_longitude_by_assay(
            dat = rgm_stage1_bAb_vac %>% 
                bind_rows(rgm_stage1_nAb_vac) %>%
                bind_rows(rgm_stage2_bAb_vac) %>%
                bind_rows(rgm_stage2_nAb_vac) %>%
                filter(!grepl("Delta", mag_cat), Group=="Vaccine") %>%
                mutate(value = mag, 
                       lower_CI=ci_l,
                       upper_CI=ci_u,
                       Ptid = Trialstage,
                       RespRate = "",
                       time = gsub(paste0(assays, collapse="|"), "", mag_cat),
                       time = factor(time, levels = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387")),
                       assay = gsub(paste0("^", times_, collapse="|"), "", mag_cat),
                       nnaive = "Non-naive",
                       lb = "",
                       lbval = 9,
                       lb2 = "",
                       lbval2 = 9) %>% 
              arrange(Trialstage, assay, time),
            x.var = "time",
            x.lb = c("D1","D22","D43","D78","D134","D202","D292","D387"),
            assays = asy,
            ylim = if (grepl("bindSpike", asy)) {c(2,6)} else if (grepl("pseudoneutid50", asy)) {c(1,5)},
            ybreaks = if (grepl("bindSpike", asy)) {2:6} else if (grepl("pseudoneutid50", asy)) {1:5},
            times = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387"),
            strip.text.x.size = 14,
            panel.text.size = 6,
            axis.text.x.size = 14,
            facet.y.var = vars(nnaive), 
            facet.x.var = vars(Trialstage),
            y.axis.lb = paste0("Geometric Mean of\n", gsub(" \\(AU/ml\\)", "", labels.assays.short[asy]))
  
        )
        
        file_name <- paste0("/", asy, "_longitudinal_gmr_nonnaive_stage1stage2.pdf")
        ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
    }
}

