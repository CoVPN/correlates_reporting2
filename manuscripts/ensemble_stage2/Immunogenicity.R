set.seed(206)
require(latex2exp)
require(vioplot)
require(lubridate)
if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
source(paste0(spath, "/code/Methods-AverageTreatmentEffect.R"))
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
ensemble.dat <- read.csv(data.path)
ensemble.dat.immuno <- read.csv(data.path)
# ensemble.dat <- read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))
# ensemble.dat.immuno <-  read.csv(paste0(spath, "/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"))

ensemble.dat <- ensemble.dat[ensemble.dat$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat$History %in% c("History a", "History b"),]
ensemble.dat <- ensemble.dat[ensemble.dat$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat$White <- with(ensemble.dat, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                  Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat$OtherRace <- with(ensemble.dat, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

ensemble.dat.immuno <- ensemble.dat.immuno[ensemble.dat.immuno$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat.immuno <- ensemble.dat.immuno[ensemble.dat.immuno$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat.immuno$History %in% c("History a", "History b") &
                                           !(is.na(ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1) & is.na(ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.1)),]
ensemble.dat.immuno$White <- with(ensemble.dat.immuno, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                                Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat.immuno$OtherRace <- with(ensemble.dat.immuno, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

# cohort x history
ensemble.dat$DelayedVaccine <- ifelse(ensemble.dat$Cohort == "Cohort 1" & ensemble.dat$History == "History a", 1, 0)
ensemble.dat$ImmediateVaccine <- ifelse(ensemble.dat$Cohort == "Cohort 1" & ensemble.dat$History == "History b", 1, 0)
ensemble.dat$DelayedHybrid <- ifelse(ensemble.dat$Cohort == "Cohort 2" & ensemble.dat$History == "History a", 1, 0)
ensemble.dat$ImmediateHybrid <- ifelse(ensemble.dat$Cohort == "Cohort 2" & ensemble.dat$History == "History b", 1, 0)

ensemble.dat$group <- character(nrow(ensemble.dat))
ensemble.dat$group <- ifelse(ensemble.dat$ImmediateVaccine == 0,
                             ensemble.dat$group, rep("Immediate Vaccine", nrow(ensemble.dat)))
ensemble.dat$group <- ifelse(ensemble.dat$DelayedVaccine == 0,
                             ensemble.dat$group, rep("Delayed Vaccine", nrow(ensemble.dat)))
ensemble.dat$group <- ifelse(ensemble.dat$ImmediateHybrid == 0,
                             ensemble.dat$group, rep("Immediate Hybrid", nrow(ensemble.dat)))
ensemble.dat$group <- ifelse(ensemble.dat$DelayedHybrid == 0,
                             ensemble.dat$group, rep("Delayed Hybrid", nrow(ensemble.dat)))

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


# weights
ensemble.dat.immuno$wt.immuno <- numeric(nrow(ensemble.dat.immuno))

# Immediate Vaccine
for(i in which(ensemble.dat.immuno$group == "Immediate Vaccine")) {
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Immediate Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Immediate Vaccine")
  }
}
# Delayed Vaccine
for(i in which(ensemble.dat.immuno$group == "Delayed Vaccine")) {
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Delayed Vaccine", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Delayed Vaccine")
  }
}
# Immediate Hybrid
for(i in which(ensemble.dat.immuno$group == "Immediate Hybrid")) {
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Immediate Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Immediate Hybrid")
  }
}
# Delayed Hybrid
for(i in which(ensemble.dat.immuno$group == "Delayed Hybrid")) {
  
  if(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjSevEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjSevEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjSevEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  } else if(ensemble.dat.immuno$AdjModEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjModEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjModEventIndPrimaryBD28_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  } else if(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc[i] == 1) {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAsymEventIndPrimaryBD28_CensorAtOutsideVacc == 1 & 
                                              ensemble.dat$AdjAsymEventIncludePrimaryBD28_CensorAtOutsideVacc == 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAsymEventIndPrimaryBD1_CensorAtOutsideVacc &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  } else {
    ensemble.dat.immuno$wt.immuno[i] <- sum(ensemble.dat$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
                                              ensemble.dat$AdjAnyEventTimePrimaryBD28_CensorAtOutsideVacc >= 1 &
                                              ensemble.dat$group == "Delayed Hybrid", na.rm = T)/
      sum(ensemble.dat.immuno$AdjAnyEventIndPrimaryBD28_CensorAtOutsideVacc == 0 &
            ensemble.dat.immuno$group == "Delayed Hybrid")
  }
}

ensemble.dat.immuno$wt.immuno <- ensemble.dat.immuno$wt.immuno/mean(ensemble.dat.immuno$wt.immuno)

## covariates
ensemble.dat.immuno$BMIDich <- ifelse(ensemble.dat.immuno$BMI %in% c(1, 2), 0, 1)
ensemble.dat.immuno$BMIDich[is.na(ensemble.dat.immuno$BMI)] <- NA

ensemble.dat.immuno$RegionUS <- ifelse(ensemble.dat.immuno$Region == 0, 1, 0)
ensemble.dat.immuno$RegionLA <- ifelse(ensemble.dat.immuno$Region == 1, 1, 0)
ensemble.dat.immuno$RegionRSA <- ifelse(ensemble.dat.immuno$Region == 2, 1, 0)

ensemble.dat.immuno$HistoryA <- ifelse(ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$CalendarDateDose2 <- as.Date(ensemble.dat.immuno$BOOSTDT)
ensemble.dat.immuno$CalendarDateDose2Early <- ifelse(ensemble.dat.immuno$CalendarDateDose2 < as.Date("2021-12-01"), 1, 0)

## covariates
ensemble.dat.immuno$BMIDich <- ifelse(ensemble.dat.immuno$BMI %in% c(1, 2), 0, 1)
ensemble.dat.immuno$BMIDich[is.na(ensemble.dat.immuno$BMI)] <- NA

ensemble.dat.immuno$RegionUS <- ifelse(ensemble.dat.immuno$Region == 0, 1, 0)
ensemble.dat.immuno$RegionLA <- ifelse(ensemble.dat.immuno$Region == 1, 1, 0)
ensemble.dat.immuno$RegionRSA <- ifelse(ensemble.dat.immuno$Region == 2, 1, 0)

ensemble.dat.immuno$HistoryA <- ifelse(ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$CalendarDateDose2 <- as.Date(ensemble.dat.immuno$BOOSTDT)
ensemble.dat.immuno$CalendarDateDose2Early <- ifelse(ensemble.dat.immuno$CalendarDateDose2 < as.Date("2021-12-01"), 1, 0)

############ Overall Marker Dist'n

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Immuno-Vioplots.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Immuno-Vioplots_", Sys.Date(), ".pdf"), width = 9, height = 7)

#BA.4.5
vioplot(ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5[ensemble.dat.immuno$group == "Immediate Vaccine"],
        ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5[ensemble.dat.immuno$group == "Delayed Vaccine"],
        ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5[ensemble.dat.immuno$group == "Immediate Hybrid"],
        ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5[ensemble.dat.immuno$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), xaxt = "n", yaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(ensemble.dat.immuno$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(ensemble.dat.immuno$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(ensemble.dat.immuno$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(ensemble.dat.immuno$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))


# B1
vioplot(ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1[ensemble.dat.immuno$group == "Immediate Vaccine"],
        ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1[ensemble.dat.immuno$group == "Delayed Vaccine"],
        ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1[ensemble.dat.immuno$group == "Immediate Hybrid"],
        ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1[ensemble.dat.immuno$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(ensemble.dat.immuno$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(ensemble.dat.immuno$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(ensemble.dat.immuno$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(ensemble.dat.immuno$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(ensemble.dat.immuno, group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                  group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                             ppr.cutoff, na.rm = TRUE),
                                digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                        group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                 na.rm = TRUE),
                                         digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                        group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                 na.rm = TRUE),
                                         digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                        group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                 na.rm = TRUE),
                                         digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(ensemble.dat.immuno,
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno,
                                                        group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                 na.rm = TRUE),
                                         digits = 1)))

dev.off()

############ Overall Marker Dist'n (USA)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Immuno-Vioplots-US.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Immuno-Vioplots-US_", Sys.Date(), ".pdf"), width = 9, height = 7)

#BA.4.5
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Vaccine"],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Vaccine"],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Hybrid"],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))


# B1
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Vaccine"],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Vaccine"],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Hybrid"],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 0)$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 0), group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 0),
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))

dev.off()

############ Overall Marker Dist'n (Latin America)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Immuno-Vioplots-LatinAmerica.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Immuno-Vioplots-LatinAmerica_", Sys.Date(), ".pdf"), width = 9, height = 7)

#BA.4.5
vioplot(subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Vaccine"],
        subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Vaccine"],
        subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Hybrid"],
        subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))


# B1
vioplot(subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Vaccine"],
        subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Vaccine"],
        subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Hybrid"],
        subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 1)$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 1), group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 1),
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))

dev.off()


############ Overall Marker Dist'n (RSA)

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Immuno-Vioplots-SouthAfrica.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Immuno-Vioplots-SouthAfrica_", Sys.Date(), ".pdf"), width = 9, height = 7)

#BA.4.5
vioplot(subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Vaccine"],
        subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Vaccine"],
        subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Hybrid"],
        subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_BA.4.5[subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_BA.4.5,
                                                   na.rm = TRUE),
                                           digits = 1)))


# B1
vioplot(subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Vaccine"],
        subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Vaccine"],
        subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Hybrid"],
        subset(ensemble.dat.immuno, Region == 2)$Year1_28dayspseudoneutid50_B.1[subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Hybrid"],
        drawRect = FALSE,
        xlab = "", ylab = "nAb titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer", border = c("black", "salmon", "blue", "gold"),
        col = c("white", "white", "white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
# axis(at = c(1, 2, 3, 4), side = 1, labels = c("Imm/Vac", "Del/Vac", "Imm/Hyb", "Del/Hyb"))
axis(at = c(1, 2, 3, 4), side = 1, labels = c("Crossover\nVaccine", "Original\nVaccine",
                                              "Crossover\nHybrid", "Original\nHybrid"),
     tick = FALSE)

points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Vaccine"), mean = 1, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Vaccine"), mean = 2, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Immediate Hybrid"), mean = 3, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")
points(x = rnorm(sum(subset(ensemble.dat.immuno, Region == 2)$group == "Delayed Hybrid"), mean = 4, sd = .1),
       y = subset(subset(ensemble.dat.immuno, Region == 2), group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Immediate Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 2, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Delayed Vaccine")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 3, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Immediate Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))
text(x = 4, y = 5, paste0("PR = ", round(100 * mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                           group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1 >
                                                      ppr.cutoff, na.rm = TRUE),
                                         digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(subset(ensemble.dat.immuno, Region == 2),
                                                          group == "Delayed Hybrid")$Year1_28dayspseudoneutid50_B.1,
                                                   na.rm = TRUE),
                                           digits = 1)))

dev.off()


############ Correlation scatterplots

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Immuno-Corplots-Jan-14-2026.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Immuno-Corplots_", Sys.Date(), ".pdf"), width = 9, height = 7)

with(subset(ensemble.dat.immuno, group == "Immediate Vaccine"),
     plot(Year1_28dayspseudoneutid50_B.1, 
          Year1_28dayspseudoneutid50_BA.4.5, 
          xlim = c(0,5), ylim = c(0,5), pch = 16, col = "gray",
          xlab = "BD-28 nAb-ID50 B.1 Titer",
          ylab = "BD-28 nAb-ID50 BA.4/5 Titer",
          main = "Crossover/Vaccine", xaxt = 'n', yaxt = 'n'))
axis(at = c(0, 1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(0, 1, 2, 3, 4), side = 1,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
with(subset(ensemble.dat.immuno, group == "Immediate Vaccine"),
     legend("topleft", paste0("Spearman correlation  = ",
                              signif(cor(Year1_28dayspseudoneutid50_B.1, Year1_28dayspseudoneutid50_BA.4.5,
                                         use = "complete.obs", method = "spearman"), digits = 3))))
smooth <- with(subset(ensemble.dat.immuno, group == "Immediate Vaccine"),
               loess(Year1_28dayspseudoneutid50_BA.4.5 ~ Year1_28dayspseudoneutid50_B.1))
lines(sort(smooth$x), smooth$fitted[order(smooth$x)], col = "salmon", lwd = 3)

with(subset(ensemble.dat.immuno, group == "Delayed Vaccine"),
     plot(Year1_28dayspseudoneutid50_B.1, 
          Year1_28dayspseudoneutid50_BA.4.5, 
          xlim = c(0,5), ylim = c(0,5), pch = 16, col = "gray",
          xlab = "BD-28 nAb-ID50 B.1 Titer",
          ylab = "BD-28 nAb-ID50 BA.4/5 Titer",
          main = "Original/Vaccine", xaxt = 'n', yaxt = 'n'))
axis(at = c(0, 1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(0, 1, 2, 3, 4), side = 1,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
with(subset(ensemble.dat.immuno, group == "Delayed Vaccine"),
     legend("topleft", paste0("Spearman correlation  = ",
                              signif(cor(Year1_28dayspseudoneutid50_B.1, Year1_28dayspseudoneutid50_BA.4.5,
                                         use = "complete.obs", method = "spearman"), digits = 3))))
smooth <- with(subset(ensemble.dat.immuno, group == "Delayed Vaccine"),
               loess(Year1_28dayspseudoneutid50_BA.4.5 ~ Year1_28dayspseudoneutid50_B.1))
lines(sort(smooth$x), smooth$fitted[order(smooth$x)], col = "salmon", lwd = 3)

with(subset(ensemble.dat.immuno, group == "Immediate Hybrid"),
     plot(Year1_28dayspseudoneutid50_B.1, 
          Year1_28dayspseudoneutid50_BA.4.5, 
          xlim = c(0,5), ylim = c(0,5), pch = 16, col = "gray",
          xlab = "BD-28 nAb-ID50 B.1 Titer",
          ylab = "BD-28 nAb-ID50 BA.4/5 Titer",
          main = "Crossover/Hybrid", xaxt = 'n', yaxt = 'n'))
axis(at = c(0, 1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(0, 1, 2, 3, 4), side = 1,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
with(subset(ensemble.dat.immuno, group == "Immediate Hybrid"),
     legend("topleft", paste0("Spearman correlation  = ",
                              signif(cor(Year1_28dayspseudoneutid50_B.1, Year1_28dayspseudoneutid50_BA.4.5,
                                         use = "complete.obs", method = "spearman"), digits = 3))))
smooth <- with(subset(ensemble.dat.immuno, group == "Immediate Hybrid"),
               loess(Year1_28dayspseudoneutid50_BA.4.5 ~ Year1_28dayspseudoneutid50_B.1))
lines(sort(smooth$x), smooth$fitted[order(smooth$x)], col = "salmon", lwd = 3)

with(subset(ensemble.dat.immuno, group == "Delayed Hybrid"),
     plot(Year1_28dayspseudoneutid50_B.1, 
          Year1_28dayspseudoneutid50_BA.4.5, 
          xlim = c(0,5), ylim = c(0,5), pch = 16, col = "gray",
          xlab = "BD-28 nAb-ID50 B.1 Titer",
          ylab = "BD-28 nAb-ID50 BA.4/5 Titer",
          main = "Original/Hybrid", xaxt = 'n', yaxt = 'n'))
axis(at = c(0, 1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(0, 1, 2, 3, 4), side = 1,
     labels = c(TeX("$10^0$"), TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
with(subset(ensemble.dat.immuno, group == "Delayed Hybrid"),
     legend("topleft", paste0("Spearman correlation  = ",
                              signif(cor(Year1_28dayspseudoneutid50_B.1, Year1_28dayspseudoneutid50_BA.4.5,
                                         use = "complete.obs", method = "spearman"), digits = 3))))
smooth <- with(subset(ensemble.dat.immuno, group == "Delayed Hybrid"),
               loess(Year1_28dayspseudoneutid50_BA.4.5 ~ Year1_28dayspseudoneutid50_B.1))
lines(sort(smooth$x), smooth$fitted[order(smooth$x)], col = "salmon", lwd = 3)

dev.off()


############ Adjusted GM (Hybrid vs Vaccine)

tab <- matrix(data = NA, nrow = 2, ncol = 4)

A <- with(ensemble.dat.immuno, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
# W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMI * Sex, BMI * (1 - Sex)))
# W <- with(ensemble.dat.immuno, cbind(RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMIDich, RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
Y.1 <- ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1
Y.45 <- ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5
weights <- ensemble.dat.immuno$wt.immuno
complete.1 <- apply(cbind(A, W, Y.1), 1, function(x){all(!is.na(x))})
complete.45 <- apply(cbind(A, W, Y.45), 1, function(x){all(!is.na(x))})
n.1 <- sum(complete.1)
n.45 <- sum(complete.45)
n <- n.1

B.1.ate <- EstimateCFMeans(A = A[complete.1], W = W[complete.1,], Y = Y.1[complete.1], weights = weights[complete.1])
B.45.ate <- EstimateCFMeans(A = A[complete.45], W = W[complete.45,], Y = Y.45[complete.45], weights = weights[complete.45])

n <- n.1
tab[1,] <- c(paste0(round(10^(B.1.ate$cf.mean.0), 2), " (", 
                   round(10^(B.1.ate$cf.mean.0 - sd(B.1.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ", ",
                   round(10^(B.1.ate$cf.mean.0 + sd(B.1.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.1.ate$cf.mean.1), 2), " (", 
                    round(10^(B.1.ate$cf.mean.1 - sd(B.1.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.1.ate$cf.mean.1 + sd(B.1.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.1.ate$ate), 2), " (", 
                    round(10^(B.1.ate$ate - sd(B.1.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.1.ate$ate + sd(B.1.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ")"),
             signif(pchisq(B.1.ate$ate^2/(var(B.1.ate$eif.ate)/n), df = 1, lower.tail = FALSE)))

n <- n.45
tab[2,] <- c(paste0(round(10^(B.45.ate$cf.mean.0), 2), " (", 
                    round(10^(B.45.ate$cf.mean.0 - sd(B.45.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.45.ate$cf.mean.0 + sd(B.45.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.45.ate$cf.mean.1), 2), " (", 
                    round(10^(B.45.ate$cf.mean.1 - sd(B.45.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.45.ate$cf.mean.1 + sd(B.45.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.45.ate$ate), 2), " (", 
                    round(10^(B.45.ate$ate - sd(B.45.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.45.ate$ate + sd(B.45.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ")"),
             signif(pchisq(B.45.ate$ate^2/(var(B.45.ate$eif.ate)/n), df = 1, lower.tail = FALSE)))

rownames(tab) <- c("B.1", "BA.45")
colnames(tab) <- c("Vaccine", "Hybrid", "Geometric Mean Ratio", "P-value")

# write.csv(tab, "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/CF-Immuno-Hybrid-vs-Vaccine.csv")
write.csv(tab, paste0(spath, "/reports/CF-Immuno-Hybrid-vs-Vaccine_", Sys.Date(), ".csv"))


############ Adjusted GM (Hybrid vs Vaccine)

tab <- matrix(data = NA, nrow = 2, ncol = 4)

A <- with(ensemble.dat.immuno, ifelse(group %in% c("Delayed Vaccine", "Delayed Hybrid"), 1, 0))
# W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMI * Sex, BMI * (1 - Sex)))
# W <- with(ensemble.dat.immuno, cbind(RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))
W <- with(ensemble.dat.immuno, cbind(Age, Sex, BMIDich, RegionUS, RegionLA, RegionRSA, CalendarDateDose2Early,
                                     with(ensemble.dat.immuno, ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))))
Y.1 <- ensemble.dat.immuno$Year1_28dayspseudoneutid50_B.1
Y.45 <- ensemble.dat.immuno$Year1_28dayspseudoneutid50_BA.4.5
weights <- ensemble.dat.immuno$wt.immuno
complete.1 <- apply(cbind(A, W, Y.1), 1, function(x){all(!is.na(x))})
complete.45 <- apply(cbind(A, W, Y.45), 1, function(x){all(!is.na(x))})
n.1 <- sum(complete.1)
n.45 <- sum(complete.45)
n <- n.1

B.1.ate <- EstimateCFMeans(A = A[complete.1], W = W[complete.1,], Y = Y.1[complete.1], weights = weights[complete.1])
B.45.ate <- EstimateCFMeans(A = A[complete.45], W = W[complete.45,], Y = Y.45[complete.45], weights = weights[complete.45])

n <- n.1
tab[1,] <- c(paste0(round(10^(B.1.ate$cf.mean.0), 2), " (", 
                    round(10^(B.1.ate$cf.mean.0 - sd(B.1.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.1.ate$cf.mean.0 + sd(B.1.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.1.ate$cf.mean.1), 2), " (", 
                    round(10^(B.1.ate$cf.mean.1 - sd(B.1.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.1.ate$cf.mean.1 + sd(B.1.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.1.ate$ate), 2), " (", 
                    round(10^(B.1.ate$ate - sd(B.1.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.1.ate$ate + sd(B.1.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ")"),
             signif(pchisq(B.1.ate$ate^2/(var(B.1.ate$eif.ate)/n), df = 1, lower.tail = FALSE)))

n <- n.45
tab[2,] <- c(paste0(round(10^(B.45.ate$cf.mean.0), 2), " (", 
                    round(10^(B.45.ate$cf.mean.0 - sd(B.45.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.45.ate$cf.mean.0 + sd(B.45.ate$eif.0)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.45.ate$cf.mean.1), 2), " (", 
                    round(10^(B.45.ate$cf.mean.1 - sd(B.45.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.45.ate$cf.mean.1 + sd(B.45.ate$eif.1)/sqrt(n) * qnorm(.975)), 2), ")"),
             paste0(round(10^(B.45.ate$ate), 2), " (", 
                    round(10^(B.45.ate$ate - sd(B.45.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ", ",
                    round(10^(B.45.ate$ate + sd(B.45.ate$eif.ate)/sqrt(n) * qnorm(.975)), 2), ")"),
             signif(pchisq(B.45.ate$ate^2/(var(B.45.ate$eif.ate)/n), df = 1, lower.tail = FALSE)))

rownames(tab) <- c("B.1", "BA.45")
colnames(tab) <- c("Crossover", "Original", "Geometric Mean Ratio", "P-value")

# write.csv(tab, "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/CF-Immuno-Original-vs-Crossover.csv")
write.csv(tab, paste0(spath, "/reports/CF-Immuno-Original-vs-Crossover_", Sys.Date(), ".csv"))
