##################################################
# obligatory to append to the top of each script #
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R")) #
##################################################

library(survey)
library(tidyverse)
library(PropCIs)
library(GGally)
library(wCorr)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
# For stratum with 1 ppt
options(survey.lonely.psu="adjust")
source("code/make_functions_graph_covail.R")
options(scipen = 999)

###################################################
#                  Parameters                     #
###################################################
# To select which tables are included in the report.
# Also to modify the headers and footers for each table.



labels.time <- c(B="Day 1",
                 Day15="Day 15",
                 Delta15overB="D15 fold-rise over D1")


labels.trt <- c("1 Dose Prototype (Moderna)", "1 Dose Beta + Omicron (Moderna)", "2 Dose Beta + Omicron (Moderna)", 
                "1 Dose Delta + Omicron (Moderna)", "1 Dose Omicron (Moderna)", "1 Dose Omicron + Prototype (Moderna)",
                "Wildtype/Prototype (Pfizer 1)", "Beta + Omicron (Pfizer 1)", "Omicron (Pfizer 1)", 
                "Beta (Pfizer 1)", "Beta + Wildtype/Prototype (Pfizer 1)","Omicron + Wildtype/Prototype (Pfizer 1)",
                "Prototype (Sanofi)", "Beta (Sanofi)", "Beta + Prototype (Sanofi)", 
                "Omicron BA.1 + Wildtype/Prototype (Pfizer 2)", "Omicron BA.4/5 + Wildtype/Prototype (Pfizer 2)" 
                )

labels.trt1 <- c("1 Dose Prototype (Moderna)", "1 Dose Beta + Omicron (Moderna)",     
                 "1 Dose Delta + Omicron (Moderna)", "1 Dose Omicron (Moderna)",         
                 "1 Dose Omicron + Prototype (Moderna)", "Prototype (Sanofi)",         
                 "Beta (Sanofi)", "Beta + Prototype (Sanofi)")

labels.trt2 <- c("Wildtype/Prototype (Pfizer 1)", "Beta + Omicron (Pfizer 1)", 
                 "Omicron (Pfizer 1)", "Beta (Pfizer 1)",  
                 "Beta + Wildtype/Prototype (Pfizer 1)", "Omicron + Wildtype/Prototype (Pfizer 1)")

if(grepl("tcell", COR)) {
  assays <- assays[substr(assays, 1, 2)=="cd"]
  } else {
  assays <- assays[substr(assays, 1, 2)!="cd"]
  }

labels.assays.short <- labels.assays.short.tabular[assays]

# redefines what is in _common.R to use shorter names
labels.assays.long <- data.frame (purrr::imap_dfc(labels.assays.short, ~ paste0(labels.assays.short[.y], ": ", labels.time)))
rownames(labels.assays.long) <- names(labels.time)

visits <- names(labels.time)#[!grepl("Delta", names(labels.time))]
assays_col <- c(config$primary_markers, config$secondary_markers, "Bcd8_FS_BA.4.5.S", "Day15cd8_FS_BA.4.5.S")

labels.assays <- expand.grid(
  time = rownames(labels.assays.long),
  marker = colnames(labels.assays.long),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    label.long = labels.assays.long[time, marker],
    label.short = sapply(labels.assays.short, as.character)[marker],
    Marker = strsplit(as.character(label.long), ": ", fixed = T)[[1]][1],
    Visit = strsplit(as.character(label.long), ": ", fixed = T)[[1]][2],
    colname = paste0(time, marker)
  )

resp.lb <- expand.grid(
  time = visits, marker = assays,
  # ind = c("Resp", "FR2", "FR4", "2lloq", "4lloq", "2llod", "4llod"), 
  ind = "Resp", 
  stringsAsFactors = F
) %>%
  mutate(Ind = case_when(
    ind == "Resp" ~ "Responder"
  )) 

labels_all <- full_join(labels.assays, resp.lb, by = c("time", "marker")) %>% 
  mutate(mag_cat = colname, resp_cat = paste0(colname, ind))

labels_all <- labels_all %>% filter(!substr(marker, nchar(marker)-1, nchar(marker)) %in% c("S1", "S2"))

###################################################
#                Clean the Data                   #
###################################################

### Table 1. Demographics 
# Output: tab_dm
# Select the covariates to be summarised.
# num_v are columns from ds_long;
# cat_v are rows of `subgroup`
tlf <-
  list(
    figs1 = list(
      fig_header = ""
    ),
    
    
    figs1supp = list(
      fig_header = ""
    ),
    
    figs2 = list(
      fig_header = ""
    ),
    
    figs3a = list(
      fig_header = ""
    ),
    
    figs3b = list(
      fig_header = ""
    ),
    
    pairs1 = list(
      fig_header = ""
    ),
    
    pairs2 = list(
      fig_header = ""
    ),
    
    pairs3 = list(
      fig_header = ""
    ),
    
    pairs4 = list(
      fig_header = ""
    ),
    
    pairs5 = list(
      fig_header = ""
    )
    
  )


for (i in names(tlf)){
  assign(i, NULL)
}

# dat.mock was made in _common.R
dat <- dat_proc


# The stratified random cohort for immunogenicity
if (is.null(dat$EthnicityUnknown)) dat$EthnicityUnknown <- 0
ds_s <- dat %>%
  mutate(
    raceC = as.character(race),
    ethnicityC = case_when(EthnicityHispanic==1 ~ "Hispanic or Latino",
                           EthnicityHispanic==0 & EthnicityNotreported==0 & 
                           EthnicityUnknown==0 ~ "Not Hispanic or Latino",
                           EthnicityNotreported==1 | 
                           EthnicityUnknown==1 ~ "Not reported and unknown "),
    RaceEthC = case_when(
      WhiteNonHispanic==1 ~ "White Non-Hispanic ",
      TRUE ~ raceC
    ),
    MinorityC = case_when(
      MinorityInd == 1 ~ "Communities of Color",
      MinorityInd == 0 ~ "White Non-Hispanic"
    ),
    AgeC = ifelse(Senior == 1, labels.age[2], labels.age[1]),
    SexC = ifelse(Sex == 1, "Female", "Male"),
    AgeSexC = paste(AgeC, SexC),
    AgeMinorC = ifelse(is.na(MinorityC), NA, paste(AgeC, MinorityC)),
    Naive = ifelse(naive==1, "Naive", "Non-naive"),
    All = "All participants"
    )

  ncol1 <- ncol(ds_s)

  ds_s <- ds_s %>% 
    mutate(Trt1=ifelse(TrtonedosemRNA==1, "One-dose mRNA booster", ""),
           Trt2=ifelse(arm %in% c(13, 14, 15), "pooled Sanofi", ""),
           Trt3=ifelse(TrtA==1, "mRNA Moderna", ""),
           Trt4=ifelse(TrtA==0, "mRNA Pfizer", ""),
           Trt5=ifelse(TrtB==1, "mRNA Prototype", ""),
           Trt6=ifelse(TrtB==1 & grepl("moderna", tolower(treatment_actual)), "mRNA Prototype Moderna", ""),
           Trt7=ifelse(TrtB==1 & grepl("pfizer", tolower(treatment_actual)), "mRNA Prototype Pfizer", ""),
           Trt8=ifelse(TrtB==0, "mRNA Omicron-Containing", ""),
           Trt9=ifelse(TrtB==0 & grepl("moderna", tolower(treatment_actual)), "mRNA Omicron-Containing Moderna", ""),
           Trt10=ifelse(TrtB==0 & grepl("pfizer", tolower(treatment_actual)), "mRNA Omicron-Containing Pfizer", ""),
           Trt11=ifelse(TrtC==1, "mRNA Bivalent", ""),
           Trt12=ifelse(TrtC==0, "mRNA Monovalent", ""),
           treatment_actual = factor(treatment_actual, levels=labels.trt)
           )
  ncol2 <- ncol(ds_s)
  Trtn <- ncol2 - ncol1
# Step2: Responders
# Post baseline visits

# assay_metadata <- assay_metadata %>% filter(panel=="tcell")
  
pos.cutoffs <- assay_metadata$pos.cutoff
names(pos.cutoffs) <- assay_metadata$assay



if(grepl("tcell", COR)){
  
  ds <- ds_s %>% 
    select(!contains("Resp", ignore.case = F)) %>%  
    select(!(contains("Day15") & contains("_resp")))
  
  names(ds) <- gsub("_resp", "Resp", names(ds))
  names(ds) <- gsub("_vacresp", "Resp", names(ds))
  
  
  
} 



plot_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x= element_text(size = 8),
        axis.text.y= element_text(size = 10),
        axis.title.x= element_text(size = 14, face="bold"),
        strip.text = element_text(size = 10, face="bold"),
        strip.background = element_rect(fill=NA,colour=NA),
        strip.placement = "outside",
        legend.position="bottom", 
        legend.text=element_text(size = 12, face="plain"),
        plot.caption=element_text(size = 11, hjust=0, face="plain"), 
        legend.box.margin=unit(c(0.2, 5.5, 0, 5.5),"inch"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 


box_assy <- c("Bcd4_IFNg.IL2_BA.4.5.S", "Day15cd4_IFNg.IL2_BA.4.5.S", "Bcd8_IFNg.IL2_BA.4.5.S", "Day15cd8_IFNg.IL2_BA.4.5.S",
              "Bcd4_IFNg.IL2_Wuhan.N", "Bcd8_IFNg.IL2_Wuhan.N", "Day15cd4_IFNg.IL2_Wuhan.N", "Day15cd8_IFNg.IL2_Wuhan.N")

ylev <- c("CD4+ T cells IFN-g/\n IL-2 BA.4.5.S (%)" ,
          "CD4+ T cells IFN-g/\n IL-2 Wuhan N (%)"  ,
          "CD8+ T cells IFN-g/\n IL-2 BA.4.5.S (%)" ,
          "CD8+ T cells IFN-g/\n IL-2 Wuhan N (%)")

chtpchs <- c(Responder=19, Negative=2)
chtcols <- c(Responder="#0AB7C9" , Negative="#8F8F8F")


# Main Figures, SubcohortInd.casedeletion

ds.i <- ds %>% dplyr::filter(!!as.name(paste0("ph1.D15", ".tcell"[grepl("tcell", COR)]))==1, SubcohortInd.casedeletion)

rpcnt_case <- get_rr(ds.i, paste0(box_assy,"Resp"), subs="All", sub.by=c("Naive", "treatment_actual"), 
                     strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 

rgm_case <- get_gm(ds.i, box_assy, subs="All", sub.by=c("Naive", "treatment_actual"), 
                   strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2, digits=3) 



ds_long <- left_join(ds %>% 
  filter(TrtonedosemRNA==1|stage==3, !!as.name(config.cor$ph2)==1, SubcohortInd.casedeletion) %>% 
  pivot_longer(cols=all_of(box_assy), values_to = "magnitude"),
  ds %>% 
    filter(TrtonedosemRNA==1|stage==3,  !!as.name(config.cor$ph2), SubcohortInd.casedeletion) %>% 
    pivot_longer(cols=all_of(paste0(grep("IFNg", box_assy, value=T), "Resp")), values_to = "response") %>% 
    mutate(name=gsub("Resp", "", name)) %>% 
    select(Ptid, name, response)
) %>% 
  left_join(rpcnt_case %>% mutate(name=gsub("Resp", "", resp_cat), rr=sprintf("%.1f%%", response*100)) %>% 
              select(Naive, name, treatment_actual, Visit, Marker, rr)) %>%
  left_join(rgm_case %>% select(Naive, name=mag_cat, treatment_actual, Visit, Marker, mag, ci_l, ci_u)) %>% 
  mutate(resp=factor(ifelse(response==0, "Negative", "Responder"), levels=c("Responder", "Negative")),
         trt_grps=ifelse(grepl("Pfizer", treatment_actual), 2, 1),
         grps=paste0(trt_grps, ifelse(grepl("CD4+", Marker), 1, 2)),
         trt_num=ifelse(trt_grps==1, match(treatment_actual, labels.trt1), match(treatment_actual, labels.trt2)),
         x_num=trt_num+ifelse(Visit=="Day 1", -1, 1)*.2 + (trt_num-1)*.05,
         group_lb=gsub("(", "\n(", treatment_actual, fixed=T),
         ytitle=gsub(" and/or", "/\n", gsub("ay ", "", sprintf("%s", Marker)), fixed=T),
         ytitle=factor(ytitle, levels=ylev),
         rslt=paste0(rr, "\n", sprintf("%.2f (%.2f, %.2f)", 10^mag, 10^ci_l, 10^ci_u))
         ) %>% 
  group_by(Naive, grps, treatment_actual) %>% 
  mutate(group_num=mean(unique(x_num))) %>% 
  ungroup()

ds_long1 <- ds_long %>% filter(!grepl("Wuhan", name)) 


figs1 <- ds_long1 %>% 
  group_split(Naive, trt_grps) %>% 
  lapply(function(x){
    
  x <- x %>% 
    filter(!is.na(response))
  
  naive <- unique(x$Naive)
  xbrks <- unique(x$x_num)
  xlbs <- paste0("D", sort(gsub("Day ", "", unique(x$Visit))))

  ggplot(x, aes(x=x_num, group=x_num, y=magnitude, color=resp))  +
    facet_grid(rows=vars(ytitle), scales = "free_y", switch = "y", labeller = labeller(groupwrap = label_wrap_gen(10))) + 
    # geom_violin(scale="width", show.legend = FALSE, color="black") +
    geom_boxplot(width=0.32, alpha = 0, outlier.shape=NA, show.legend = FALSE, color="black") +
    geom_jitter(aes(color=resp, shape=resp), width = 0.15, height = 0, alpha=0.75, size =1.2, show.legend = TRUE) +
    geom_text(aes(label = rslt, y = Inf), vjust = 1.5, size=2.4, color="black",  check_overlap = T) +
    geom_text(aes(x=group_num, y = -Inf, label=group_lb), size=2.8, vjust=1.8, color="black", check_overlap = T) +
    labs(x="Vaccine Arm", title=sprintf("%s Cohort", naive), y="") +
    coord_cartesian(clip="off") +
    scale_shape_manual(name="", values=chtpchs) +
    scale_color_manual(name="", values=chtcols) +
    scale_x_continuous(breaks = xbrks, 
                       labels=rep(xlbs, length(xbrks)/2),
                       expand = expansion(add = c(.1, .1))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, .25)), 
                       breaks = log10(c(0.002, .01, .05, .1, .2,  .5,  1, 2, 5, 10, 20)),
                       labels = paste0(c(0.002, .01, .05, .1, .2, .5, 1, 2, 5, 10, 20), "%"))+
    plot_theme +
    theme(axis.title.x= element_text(vjust = -6))
  })


# Supp Figures 
# S2-3, SubcohortInd.casedeletion


ds_long2 <- ds_long %>% 
  filter(grepl("Wuhan",name)) %>% 
  mutate(rslt=paste0(rr, "\n", sprintf("%.3f (%.3f, %.3f)", 10^mag, 10^ci_l, 10^ci_u)),
         x_num=trt_num+ifelse(Visit=="Day 1", -1, 1)*.22 + (trt_num-1)*.02) %>% 
  group_by(Naive, grps, treatment_actual) %>% 
  mutate(group_num=mean(unique(x_num))) %>% 
  ungroup()


figs1supp <- 
  ds_long2 %>% 
  group_split(Naive, trt_grps) %>% 
  lapply(function(x){
    x <- x %>% 
      filter(!is.na(response))
    
    naive <- unique(x$Naive)
    xbrks <- unique(x$x_num)
    xlbs <- paste0("D", sort(gsub("Day ", "", unique(x$Visit))))
    
    ggplot(x, aes(x=x_num, group=x_num, y=magnitude, color=resp))  +
      facet_grid(rows=vars(ytitle), scales = "free_y", switch = "y", labeller = labeller(groupwrap = label_wrap_gen(10))) + 
      # geom_violin(scale="width", show.legend = FALSE, color="black") +
      geom_boxplot(width=0.32, alpha = 0, outlier.shape=NA, show.legend = FALSE, color="black") +
      geom_jitter(aes(color=resp, shape=resp), width = 0.15, height = 0, alpha=0.75, size =1.2, show.legend = TRUE) +
      geom_text(aes(label = rslt, y = Inf), vjust = 1.5, size=2.15, color="black",  check_overlap = T) +
      geom_text(aes(x=group_num, y = -Inf, label=group_lb), size=2.8, vjust=1.8, color="black", check_overlap = T) +
      labs(x="Vaccine Arm", title=sprintf("%s Cohort", naive), y="") +
      coord_cartesian(clip="off") +
      scale_shape_manual(name="", values=chtpchs) +
      scale_color_manual(name="", values=chtcols) +
      scale_x_continuous(breaks = xbrks, 
                         labels=rep(xlbs, length(xbrks)/2),
                         expand = expansion(add = c(.1, .1))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, .25)), 
                         breaks = log10(c(0.002, .01, .05, .1, .2,  .5,  1, 2, 5, 10, 20)),
                         labels = paste0(c( 0.002, .01, .05, .1, .2, .5, 1, 2, 5, 10, 20), "%"))+
      plot_theme +
      theme(axis.title.x= element_text(vjust = -6))
  })


ylev_fs <- c("CD4+ T cells FS BA.4.5.S", "CD8+ T cells FS BA.4.5.S" )

rgm_fs <- get_gm(ds.i, grep("BA.4.5.S", FS_vars, fixed=T, value=T), subs="All", 
                 sub.by=c("Naive", "treatment_actual"), 
                 strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2, digits=2) 


ds_fs <- ds %>% 
  filter(TrtonedosemRNA==1|stage==3, !!as.name(config.cor$ph2)==1, SubcohortInd.casedeletion) %>% 
  pivot_longer(cols=all_of(grep("BA.4.5.S", FS_vars, fixed=T, value=T)), values_to = "fs") %>% 
  left_join(labels_all %>% rename(name=colname)) %>% 
  left_join(rgm_fs %>% select(Naive, name=mag_cat, treatment_actual, Visit, Marker, mag, ci_l, ci_u)) %>% 
  mutate(ytitle=gsub("functionality score", "FS", gsub("ay ", "", sprintf("%s", Marker)), fixed=T),
         ytitle=factor(ytitle, levels=ylev_fs),
         trt_grps=ifelse(grepl("Pfizer", treatment_actual), 2, 1),
         grps=paste0(trt_grps, ifelse(grepl("CD4+", Marker), 1, 2)),
         trt_num=ifelse(trt_grps==1, match(treatment_actual, labels.trt1), match(treatment_actual, labels.trt2)),
         x_num=trt_num+ifelse(Visit=="Day 1", -1, 1)*.22 + (trt_num-1)*.05,
         group_lb=gsub("(", "\n(", treatment_actual, fixed=T),
         rslt=sprintf("%.3f\n(%.3f, %.3f)", mag, ci_l, ci_u)
  ) %>% 
  group_by(Naive, grps, treatment_actual) %>% 
  mutate(group_num=mean(unique(x_num))) %>% 
  ungroup()

ds_fs <- ds_fs %>%
  mutate(cohort_event = case_when(COVIDIndD22toD91==1 ~ "Booster-Proximal Cases", #"D22-D91 Cases",
                                  COVIDIndD92toD181==1 ~ "Booster-Distal Cases", #"D92-D181 Cases",
                                  COVIDIndD22toD181==0 ~ "Non-Cases")) %>%
  rbind(ds_fs %>%
          filter(COVIDIndD22toD181==1) %>% 
          mutate(cohort_event="Cases (Proximal + Distal)")) %>% #"D22-D181 Cases"
  mutate(cohort_event=factor(cohort_event,levels = c("Non-Cases","Booster-Proximal Cases","Booster-Distal Cases","Cases (Proximal + Distal)")))


figs2 <- 
  ds_fs %>% 
  group_split(Naive, trt_grps) %>% 
  lapply(function(x){
    x <- x %>% 
      filter(!is.na(fs))
    naive <- unique(x$Naive)
    xbrks <- unique(x$x_num)
    xlbs <- paste0("D", sort(gsub("Day ", "", unique(x$Visit))))
    
    
    ggplot(x, aes(x=x_num, group=x_num, y=fs))  +
      facet_grid(rows=vars(ytitle), scales = "free_y", switch = "y", labeller = labeller(groupwrap = label_wrap_gen(10))) + 
      # geom_violin(scale="width", color="#0AB7C9", show.legend = FALSE) +
      geom_boxplot(width=0.32, color="#0AB7C9", alpha = 0, outlier.shape=NA, show.legend = FALSE) +
      geom_jitter(width = 0.15, color="#0AB7C9", height = 0, alpha=0.5, size =1.2, show.legend = FALSE) +
      geom_text(aes(label = rslt, y = Inf), vjust = 1.5, size=2.2, color="black",  check_overlap = T) +
      geom_text(aes(x=group_num, y = -Inf, label=group_lb), size=2.8, vjust=1.8, color="black", check_overlap = T) +
      labs(x="", title=sprintf("Distributions of T-cell functionality score markers for each one-dose vaccine booster arm\n%s Cohort", naive), y="") +
      coord_cartesian(clip="off") +
      scale_x_continuous(breaks = xbrks, 
                         labels=rep(xlbs, length(xbrks)/2),
                         expand = expansion(add = c(.1, .1))) +
      scale_y_continuous(expand = expansion(mult = c(0.06, .2))) +
      plot_theme 
  })

# Supp 4-5

more_assy <- c("Bcd4_IFNg.IL2_BA.4.5.S", "Bcd4_FS_BA.4.5.S", "Day15cd4_FS_BA.4.5.S", "Day15cd4_IFNg.IL2_BA.4.5.S",
               "Bcd8_IFNg.IL2_BA.4.5.S", "Bcd8_FS_BA.4.5.S", "Day15cd8_FS_BA.4.5.S", "Day15cd8_IFNg.IL2_BA.4.5.S",
               "Bcd4_IFNg.IL2_Wuhan.N", "Bcd8_IFNg.IL2_Wuhan.N")

ylev_more <- c("D1 CD4+ T cells IFN-g/\n IL-2 BA.4.5.S (%)" ,
               "D1 CD4+ T cells FS BA.4.5.S", 
               "D15 CD4+ T cells IFN-g/\n IL-2 BA.4.5.S (%)",
               "D15 CD4+ T cells FS BA.4.5.S",
               "D1 CD8+ T cells IFN-g/\n IL-2 BA.4.5.S (%)" ,
               "D1 CD8+ T cells FS BA.4.5.S", 
               "D15 CD8+ T cells IFN-g/\n IL-2 BA.4.5.S (%)",
               "D15 CD8+ T cells FS BA.4.5.S",
               "D1 CD4+ T cells IFN-g/\n IL-2 Wuhan N (%)",
               "D1 CD8+ T cells IFN-g/\n IL-2 Wuhan N (%)")

ds_long_more <- left_join(ds %>% 
                       filter(TrtonedosemRNA==1, ph2.D15.tcell) %>% 
                       pivot_longer(cols=all_of(more_assy), values_to = "magnitude"),
                     ds %>% 
                       filter(TrtonedosemRNA==1, ph2.D15.tcell) %>% 
                       pivot_longer(cols=all_of(paste0(grep("IFNg", box_assy, value=T), "Resp")), values_to = "response") %>% 
                       mutate(name=gsub("Resp", "", name)) %>% 
                       select(Ptid, name, response)) 

ds_long_case1 <- ds_long_more %>%
  filter(TrtonedosemRNA==1) %>% 
  mutate(cohort_event = case_when(COVIDIndD22toD91==1 ~ "Booster-Proximal Cases", #"D22-D91 Cases",
                                  COVIDIndD92toD181==1 ~ "Booster-Distal Cases", #"D92-D181 Cases",
                                  COVIDIndD22toD181==0 ~ "Non-Cases"))
ds_long_case2 <- ds_long_more %>%
  filter(TrtonedosemRNA==1) %>% 
  mutate(cohort_event = case_when(COVIDIndD22toD181==1 ~ "Cases (Proximal + Distal)", #"D22-D181 Cases",
                                  COVIDIndD22toD181==0 ~ "Non-Cases"))
  

ds.i.long1 <- ds_long_case1 %>% dplyr::filter(!!as.name(config.cor$ph1)==1)
ds.i.long2 <- ds_long_case2 %>% dplyr::filter(!!as.name(config.cor$ph1)==1)

rpcnt_case1 <- get_rr(ds.i.long1,
                      paste0(box_assy,"Resp"), 
                      subs="cohort_event", 
                      sub.by=c("Naive"), 
                      strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 


rpcnt_case2 <- get_rr(ds.i.long2,
                      paste0(box_assy,"Resp"), 
                      subs="cohort_event", 
                      sub.by=c("Naive"), 
                      strata=config.cor$WtStratum, weights=config.cor$wt, subset=config.cor$ph2) 

rpcnt_case_long <- full_join(rpcnt_case1, rpcnt_case2)


ds_long_case <- bind_rows(ds_long_case1, 
                          ds_long_case2 %>% filter(cohort_event=="Cases (Proximal + Distal)")) %>% 
  left_join(rpcnt_case_long %>% 
              mutate(name=gsub("Resp", "", resp_cat), rr=sprintf("%.1f%%", response*100)) %>%
              select(Naive, name, cohort_event=Group, rr)) %>%  
  left_join(labels_all %>% select(name=colname, Visit, Marker)) %>% 
  mutate(resp=factor(ifelse(grepl("FS", name), cohort_event, ifelse(response==0, "Negative", cohort_event))),
         ytitle=gsub("functionality score", "FS", gsub(" and/or", "/\n", gsub("ay ", "", sprintf("%s %s", Visit, Marker)), fixed=T)),
         ytitle=factor(ytitle, levels=ylev_more),
         tcell=ifelse(grepl("cd4_", name), "CD4", "CD8"),
         cohort_event=factor(cohort_event,levels = c("Non-Cases","Booster-Proximal Cases","Booster-Distal Cases","Cases (Proximal + Distal)"))) %>% 
  filter(Naive=="Non-naive"|(Naive=="Naive" & !name%in%c("Bcd4_IFNg.IL2_Wuhan.N", "Bcd8_IFNg.IL2_Wuhan.N")))


casecols <- c(Negative="#8F8F8F", `Non-Cases`="#D92321", `Booster-Proximal Cases`="#0AB7C9", 
              `Booster-Distal Cases`="#0AB7C9", `Cases (Proximal + Distal)`="#0AB7C9")

casepchs <- c(Negative=2, `Non-Cases`=19, `Booster-Proximal Cases`=19, 
              `Booster-Distal Cases`=19, `Cases (Proximal + Distal)`=19)

figs3a <- 
  ds_long_case %>% 
  filter(!grepl("FS", name)) %>% 
  group_split(Naive) %>% 
  lapply(function(x){
    x <- x %>% 
      filter(!is.na(resp))
    naive <- unique(x$Naive)
    nrow <- n_distinct(x$ytitle)
    
    ggplot(x, aes(x=cohort_event, y=magnitude, color=cohort_event))  +
      facet_grid(rows=vars(ytitle), scales = "free_y", switch = "y", labeller = labeller(groupwrap = label_wrap_gen(10))) + 
      geom_violin(scale="width", show.legend = FALSE) +
      geom_boxplot(width=0.25, alpha = 0, outlier.shape=NA, show.legend = FALSE) +
      geom_jitter(aes(color=resp, shape=resp), width = 0.1, height = 0, alpha=0.5, size =1.2, show.legend = FALSE) +
      geom_text(aes(label = rr, y = Inf), vjust = 1.5, size=3, color="black", position = position_dodge(width = .8), check_overlap = T) +
      labs(title=sprintf("Distributions of T-cell responses among the one-dose mRNA group\n%s Cohort", naive), y="", x="") +
      scale_shape_manual(name="", values=casepchs) +
      scale_color_manual(name="", values=casecols) +
      scale_x_discrete(labels = scales::label_wrap(11)) +
      scale_y_continuous(expand = expansion(mult = c(.06, .25)), 
                         breaks = log10(c(0.002, .01, .1, .2,  .5,  1, 2, 5, 10, 20)),
                         labels = paste0(c(0.002, .01, .1, .2, .5, 1, 2, 5, 10, 20), "%")) +
      plot_theme +
      theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=ifelse(nrow>4, 8, 10)))
  })


figs3b <- 
  ds_long_case %>% 
  filter(grepl("FS", name)) %>% 
  group_split(Naive) %>% 
  lapply(function(x){
    x <- x %>% 
      filter(!is.na(resp))
    
    naive <- unique(x$Naive)
    ggplot(x, aes(x=cohort_event, y=magnitude, color=cohort_event))  +
      facet_grid(rows=vars(ytitle), scales = "free_y", switch = "y", labeller = labeller(groupwrap = label_wrap_gen(10))) + 
      geom_violin(scale="width", show.legend = FALSE) +
      geom_boxplot(width=0.25, alpha = 0, outlier.shape=NA, show.legend = FALSE) +
      geom_jitter(aes(color=resp, shape=resp), width = 0.1, height = 0, alpha=0.5, size =1.2, show.legend = FALSE) +
      # geom_text(aes(label = rr, y = Inf), vjust = 1.5, size=3, color="black", position = position_dodge(width = .8), check_overlap = T) +
      labs(title=sprintf("Distributions of T-cell responses among the one-dose mRNA group\n%s Cohort", naive), y="", x="") +
      scale_shape_manual(name="", values=casepchs) +
      scale_color_manual(name="", values=casecols) +
      scale_x_discrete(labels = scales::label_wrap(11)) +
      scale_y_continuous(expand = expansion(mult = c(.06, .25))) +
      plot_theme +
      theme(axis.text.x = element_text(size=10))
  })


# Supp 6

  plot_dat <- ds %>% dplyr::filter(!!as.name(paste0("ph1.D15", ".tcell"[grepl("tcell", COR)]))==1) 
  strata <- config.cor$WtStratum
  weight <- config.cor$wt
  plot_title <- "Correlations of T-cell and antibody marker in the Per-Protocol Immunogenicity Analysis Set"
  
  corcols <- c("Day15pseudoneutid50_BA.4.BA.5", 
               "Bcd4_IFNg.IL2_BA.4.5.S", "Day15cd4_IFNg.IL2_BA.4.5.S", 
               "Bcd8_IFNg.IL2_BA.4.5.S", "Day15cd8_IFNg.IL2_BA.4.5.S",
               "Bcd4_IFNg.IL2_Wuhan.N", "Day15cd4_IFNg.IL2_Wuhan.N", 
               "Bcd8_IFNg.IL2_Wuhan.N", "Day15cd8_IFNg.IL2_Wuhan.N",
               "Bcd4_IFNg.IL2_COV2.CON.S", "Day15cd4_IFNg.IL2_COV2.CON.S", 
               "Bcd8_IFNg.IL2_COV2.CON.S", "Day15cd8_IFNg.IL2_COV2.CON.S")
  
  corlabels <- c("D15 nAb-ID50 BA.4/5", 
                 "D1 CD4+ IFN-γ/IL-2 Spike BA.4/5", "D15 CD4+ IFN-γ/IL-2 Spike BA.4/5",
                 "D1 CD8+ IFN-γ/IL-2 Spike BA.4/5", "D15 CD8+ IFN-γ/IL-2 Spike BA.4/5",
                 "D1 CD4+ IFN-γ/IL-2 N Index", "D15 CD4+ IFN-γ/IL-2 N Index",
                 "D1 CD8+ IFN-γ/IL-2 N Index", "D15 CD8+ IFN-γ/IL-2 N Index",
                 "D1 CD4+ IFN-γ/IL-2 Spike Index", "D15 CD4+ IFN-γ/IL-2 Spike Index",
                 "D1 CD8+ IFN-γ/IL-2 Spike Index", "D15 CD8+ IFN-γ/IL-2 Spike Index")

  names(corlabels) <- corcols


pairs1 <- pairplots(data=plot_dat, plot_title=plot_title, strata=strata, weight=weight, 
					corcols=c("Day15pseudoneutid50_BA.4.BA.5", 
							  "Bcd4_IFNg.IL2_BA.4.5.S", "Day15cd4_IFNg.IL2_BA.4.5.S", "Bcd4_IFNg.IL2_Wuhan.N",
							  "Bcd8_IFNg.IL2_BA.4.5.S", "Day15cd8_IFNg.IL2_BA.4.5.S"))

pairs2 <- pairplots(data=plot_dat, plot_title="Correlations of CD4+ T-cell marker in the Per-Protocol Immunogenicity Analysis Set",
                    strata=strata, weight=weight, 
                    corcols=c("Bcd4_IFNg.IL2_BA.4.5.S", "Day15cd4_IFNg.IL2_BA.4.5.S", 
                              "Bcd4_IFNg.IL2_Wuhan.N", "Day15cd4_IFNg.IL2_Wuhan.N"))
  
pairs3 <- pairplots(data=plot_dat, plot_title="Correlations of CD8+ T-cell marker in the Per-Protocol Immunogenicity Analysis Set", 
                    strata=strata, weight=weight, 
                    corcols=c("Bcd8_IFNg.IL2_BA.4.5.S", "Day15cd8_IFNg.IL2_BA.4.5.S", 
                              "Bcd8_IFNg.IL2_Wuhan.N", "Day15cd8_IFNg.IL2_Wuhan.N"))



pairs4 <- pairplots(data=plot_dat, plot_title="Correlations of CD4+ T-cell marker in the Per-Protocol Immunogenicity Analysis Set", 
                    strata=strata, weight=weight, 
                    corcols=c("Bcd4_IFNg.IL2_BA.4.5.S", "Day15cd4_IFNg.IL2_BA.4.5.S", 
                              "Bcd4_IFNg.IL2_COV2.CON.S", "Day15cd4_IFNg.IL2_COV2.CON.S"))


pairs5 <- pairplots(data=plot_dat, plot_title="Correlations of CD8+ T-cell marker in the Per-Protocol Immunogenicity Analysis Set", 
                    strata=strata, weight=weight, 
                    corcols=c("Bcd8_IFNg.IL2_BA.4.5.S", "Day15cd8_IFNg.IL2_BA.4.5.S", 
                              "Bcd8_IFNg.IL2_COV2.CON.S", "Day15cd8_IFNg.IL2_COV2.CON.S"))


# path for tables
save.results.to <- here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to <- paste0(here::here("output"), "/", attr(config,"config"))

if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

save(list = c("tlf", names(tlf)), 
     file = file.path(save.results.to, sprintf("Figure%s.Rdata", ifelse(exists("COR"), gsub("D15", "D22", COR), ""))))

