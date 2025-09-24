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
    
}

    
  for (i in c("bAb", "nAb")){
    for (j in 1:2){
    ylim <- if(i=="bAb"){c(2,6)}else{c(1,5)}
    assys <- if(i=="bAb"){c("bindSpike", "bindSpike_beta", "bindSpike_omicron")
    }else{
        c("pseudoneutid50", "pseudoneutid50_B.1.351", "pseudoneutid50_BA.1")}
    asy <- ifelse(i=="bAb", "Binding Antibody", "Psv Neutralization Antibody")
    y.lb.scale <- "log"
    
    dat_plot <- get(sprintf("rgm_stage%s_%s_vac", j,i)) %>% 
                    rowwise() %>%
                    mutate(value = mag,
                           lower_CI=ci_l,
                           upper_CI=ci_u,
                           Ptid = Trialstage,
                           RespRate = "",
                           time = gsub(paste0(assays, collapse="|"), "", mag_cat),
                           time = factor(time, levels = c("B","Day22","Day43","Day78","Day134","Day202","Day292","Day387")),
                           marker = gsub(time, "", mag_cat),
                           assay = gsub(paste0("^", times_, collapse="|"), "", mag_cat),
                           nnaive = "Non-naive",
                           lb = "",
                           lbval = 9,
                           lb2 = "",
                           lbval2 = 9) %>% 
      filter(assay %in% assys & 
               time %in% c("B","Day22","Day43","Day78","Day134","Day202")) %>%
      left_join(assay_metadata, by="assay") %>%
      mutate(panel = ifelse(grepl("pseudo", assay), "nAb ID50", ifelse(grepl("bindSpike", assay), "Binding IgG Spike", ""))) %>%
      mutate(assay_label_short = gsub("PsV Neutralization to |PsV Neutralization |Binding Antibody to Spike ", "", assay_label)) %>%
      ungroup()
    

      f_2 <- dat_plot %>%    
        mutate(assay_label=factor(assay_label, levels=c("Binding Antibody to Spike D614",
                                                        "Binding Antibody to Spike Beta",
                                                        "Binding Antibody to Spike BA.1", 
                                                        "PsV Neutralization to D614G",                    
                                                        "PsV Neutralization to Beta",                     
                                                        "PsV Neutralization to BA.1"      ))) %>% 
      ggplot(aes(x = time, y = value)) +
      facet_grid(col = vars(assay_label)) +
      geom_line(aes(group = Group, linetype=Group, color = Group)) +
      scale_linetype_manual(values= c("Vaccine"=1, "Placebo"=2)) + 
      geom_point(aes(color=Group), size = 3, alpha = 0.6, show.legend = TRUE) +
      geom_text(aes(x = time, label = !!sym("RespRate"), y = ylim[2]*0.9), color = "black", size = 16, check_overlap = TRUE) +
      geom_hline(aes(yintercept = ifelse(RespRate!="",lbval,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
      geom_text(aes(label = ifelse(RespRate!="",lb,""), x = 0.4, y = lbval), hjust = 0, color = "black", size = 16, check_overlap = TRUE, na.rm = TRUE) + 
      geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI, linetype=Group, color=Group), width=.2) +
      geom_hline(aes(yintercept = ifelse(RespRate!="",lbval2,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
      geom_text(aes(label = ifelse(RespRate!="",lb2,""), x = 0.4, y = lbval2), hjust = 0, color = "black", size = 16, check_overlap = TRUE, na.rm = TRUE) + 
      scale_color_manual(values=c("Vaccine"="black", "Placebo"="grey55")) +
      scale_x_discrete(labels = c("D1","D22","D43","D78","D134","D202","D292","D387"), drop=TRUE) +
      scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2]), labels = ifelse(y.lb.scale == "log", scales::math_format(10^.x))) +
      labs(x = "Time", y = sprintf("Geometric Mean of %s", asy), 
      title =sprintf("Geometric Mean of %s longitudinal plots across timepoints\nStage %s", asy, j), shape="Treatment Arm") +
      theme_bw(base_size = 16) +
      theme(
            strip.background = element_rect(fill=NA,colour=NA),
            strip.placement = "outside",
            legend.position = "bottom",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
    
              file_name <- sprintf("/%s_longitudinal_gmr_nonnaive_stage%s.pdf", asy, j)
              ggsave(plot = f_2, filename = paste0(save.results.to, file_name), width = 16, height = 11)
      
    }
  }  
    