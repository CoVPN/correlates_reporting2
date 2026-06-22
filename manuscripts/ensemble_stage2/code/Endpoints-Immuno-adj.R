require(lubridate)
require(survminer)
require(survival)
set.seed(206)

if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
# ensemble.dat <- read.csv(data.path) # read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))
ensemble.dat.immuno <- read.csv(data.path) #read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))
ensemble.dat.immuno <- ensemble.dat.immuno[ensemble.dat.immuno$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat.immuno <- ensemble.dat.immuno[ensemble.dat.immuno$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat.immuno$History %in% c("History a", "History b") &
                                             !(is.na(ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1) & is.na(ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.1)),]
ensemble.dat.immuno$White <- with(ensemble.dat.immuno, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                  Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat.immuno$OtherRace <- with(ensemble.dat.immuno, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

# cohort x history
ensemble.dat.immuno$DelayedVaccine <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 1" & ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$ImmediateVaccine <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 1" & ensemble.dat.immuno$History == "History b", 1, 0)
ensemble.dat.immuno$DelayedHybrid <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 2" & ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$ImmediateHybrid <- ifelse(ensemble.dat.immuno$Cohort == "Cohort 2" & ensemble.dat.immuno$History == "History b", 1, 0)

ensemble.dat.immuno$group <- character(nrow(ensemble.dat.immuno))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$ImmediateVaccine == 0,
                                    ensemble.dat.immuno$group, rep("Immediate Vaccine", nrow(ensemble.dat.immuno)))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$DelayedVaccine == 0,
                                    ensemble.dat.immuno$group, rep("Delayed Vaccine", nrow(ensemble.dat.immuno)))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$ImmediateHybrid == 0,
                                    ensemble.dat.immuno$group, rep("Immediate Hybrid", nrow(ensemble.dat.immuno)))
ensemble.dat.immuno$group <- ifelse(ensemble.dat.immuno$DelayedHybrid == 0,
                                    ensemble.dat.immuno$group, rep("Delayed Hybrid", nrow(ensemble.dat.immuno)))


### Event variables for correlates analyses

## day 28
D28.gap.immuno <- as.numeric(as.Date(ensemble.dat.immuno$SVSTDTC_YEAR.1_28.DAYS)) - 
  as.numeric(as.Date(ensemble.dat.immuno$BOOSTDT))
D28.gap <- as.numeric(as.Date(ensemble.dat$SVSTDTC_YEAR.1_28.DAYS)) - 
  as.numeric(as.Date(ensemble.dat$BOOSTDT))


# correlates cohort
ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno < 1,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc)


# per-protocol
ensemble.dat$AdjModtoSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap < 1,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap)
ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventIncludePrimaryBD28_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc == 0, 0, 1)



## day 35

# correlates cohort
ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc)

ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno <= 7,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap.immuno - 7)
ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc)


# per-protocol
ensemble.dat$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModtoSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjModtoSevEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjModEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjModEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjSevEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjSevEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjAsymEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)

ensemble.dat$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc == 1 & 
           ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAnyEventIndPrimaryBD1_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap <= 7,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD1_CensorAtOutsideVacc - D28.gap - 7)
ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <- 
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc < 0,
         0, ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc)
ensemble.dat$AdjAnyEventIncludePrimaryBD35_CensorAtOutsideVacc <-
  ifelse(ensemble.dat$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc == 0, 0, 1)


########################################
### Event counts and event rate by group
########################################

Tab <- matrix(nrow = 6, ncol = 5)
colnames(Tab) <- c("Crossover/Vaccine", "Crossover/Hybrid", "Original/Vaccine", "Original/Hybrid", "Overall")
rownames(Tab) <- c("N", "Moderate to Severe-critical COVID (Total)",
                   "Moderate COVID (Total)",
                   "Severe-critical COVID (Total)",
                   "Asymptomatic Infection (Total)",
                   "Infection, Any Symptomatology (Total)")
# N
Tab[1,] <- c(sum(ensemble.dat.immuno$ImmediateVaccine), sum(ensemble.dat.immuno$ImmediateHybrid),
              sum(ensemble.dat.immuno$DelayedVaccine), sum(ensemble.dat.immuno$DelayedHybrid), nrow(ensemble.dat.immuno))

# Moderate to Severe-critical COVID
event.total <- with(ensemble.dat.immuno, 
                    c(sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[2,] <- event.total
# annual.incidence <- with(ensemble.dat.immuno,
#                        c(sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                            sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                          sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                            sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                          sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                            sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                          sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                            sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                          sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                            sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[2,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")


# Moderate COVID
ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc_ME <- 
  ifelse(ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
         ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1,
         0, ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc)
event.total <- with(ensemble.dat.immuno, 
                    c(sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc_ME[ImmediateVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc_ME[ImmediateHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc_ME[DelayedVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc_ME[DelayedHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc_ME)))
Tab[3,] <- event.total

# annual.incidence <- with(ensemble.dat.immuno,
#                          c(sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[3,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Severe-critical COVID
event.total <- with(ensemble.dat.immuno, 
                    c(sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[4,] <- event.total

# annual.incidence <- with(ensemble.dat.immuno,
#                          c(sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[4,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Asymp
#set event indicator to 0 when asym and mod-severe covid occur concurrently
ensemble.dat.immuno$AdjAsymOtherEventIndPrimaryBD35_CensorAtOutsideVacc_ME <- 
  ifelse(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
         ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1,
         0, ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc)
event.total <- with(ensemble.dat.immuno, 
                    c(sum(AdjAsymOtherEventIndPrimaryBD35_CensorAtOutsideVacc_ME[ImmediateVaccine == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD35_CensorAtOutsideVacc_ME[ImmediateHybrid == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD35_CensorAtOutsideVacc_ME[DelayedVaccine == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD35_CensorAtOutsideVacc_ME[DelayedHybrid == 1]),
                      sum(AdjAsymOtherEventIndPrimaryBD35_CensorAtOutsideVacc_ME)))
Tab[5,] <- event.total


# annual.incidence <- with(ensemble.dat.immuno,
#                          c(sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[5,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Infection, Any Symptomatology
event.total <- with(ensemble.dat.immuno, 
                    c(sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[6,] <- event.total

# annual.incidence <- with(ensemble.dat.immuno,
#                          c(sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[6,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

write.csv(Tab, paste0(spath, "/reports/AdjEndpointsImmuno_", Sys.Date(),".csv"))  


########################################
### Event counts and event rate by group
########################################

Tab <- matrix(nrow = 6, ncol = 5)
colnames(Tab) <- c("Immediate Vacc", "Immedate Hybrid", "Delayed Vacc", "Delayed Hybrid", "Overall")
rownames(Tab) <- c("N", "Moderate to Severe-critical COVID (Total)",
                   "Moderate COVID (Total)",
                   "Severe-critical COVID (Total)",
                   "Asymptomatic Infection (Total)",
                   "Infection, Any Symptomatology (Total)")
# N
Tab[1,] <- c(sum(ensemble.dat.immuno$ImmediateVaccine), sum(ensemble.dat.immuno$ImmediateHybrid),
             sum(ensemble.dat.immuno$DelayedVaccine), sum(ensemble.dat.immuno$DelayedHybrid), nrow(ensemble.dat.immuno))

# Moderate to Severe-critical COVID
event.total <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
                                   ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[2,] <- event.total


# annual.incidence <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
#                                         ensemble.dat.immuno$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 90),],
#                          c(sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[2,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")


# Moderate COVID
event.total <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
                                   ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[3,] <- event.total


# annual.incidence <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjModEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
#                                         ensemble.dat.immuno$AdjModEventTimePrimaryBD35_CensorAtOutsideVacc > 90),],
#                          c(sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjModEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjModEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[3,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Severe-critical COVID
event.total <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
                                   ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[4,] <- event.total


# annual.incidence <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
#                                         ensemble.dat.immuno$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > 90),],
#                          c(sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[4,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Asymp
event.total <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
                                   ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[5,] <- event.total


# annual.incidence <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
#                                         ensemble.dat.immuno$AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc > 90),],
#                          c(sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjAsymEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjAsymEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[5,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# Infection, Any Symptomatology
event.total <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
                                   ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 90),], 
                    c(sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
                      sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc)))
Tab[6,] <- event.total

# annual.incidence <- with(ensemble.dat.immuno[!(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
#                                         ensemble.dat.immuno$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 90),],
#                          c(sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateVaccine == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[ImmediateHybrid == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedVaccine == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1])/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[DelayedHybrid == 1]),
#                            sum(AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc)/
#                              sum(AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc))) * 365
# Tab[6,] <- paste0(round(event.total, digits = 1), " (", round(100 * annual.incidence, digits = 1), "%)")

# write.csv(Tab, "/Volumes/ahudson/Ensemble-stage2-correlates/SummaryTables/AdjEndpointsImmuno90Days.csv")
write.csv(Tab, paste0(spath, "/reports/AdjEndpointsImmuno90Days_", Sys.Date(),".csv"))  
