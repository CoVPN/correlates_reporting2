set.seed(206)
require(latex2exp)
require(vioplot)
require(lubridate)
if(!exists("spath")) {
  spath <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2"
}
source(paste0(spath, "/code/Methods-AverageTreatmentEffect.R"))
source(paste0(spath, "/code/Methods-NPCorrelates.R"))
if(!exists("data.path")) {
  data.path <- "/Volumes/trials/covpn/p3003/analysis/correlates/stage2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv"
}
ensemble.dat <- read.csv(data.path)
ensemble.dat.immuno <- read.csv(data.path)

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

ensemble.dat.immuno$HistoryA <- ifelse(ensemble.dat.immuno$History == "History a", 1, 0)
ensemble.dat.immuno$CalendarDateDose2 <- as.Date(ensemble.dat.immuno$BOOSTDT)
ensemble.dat.immuno$CalendarDateDose2Early <- ifelse(ensemble.dat.immuno$CalendarDateDose2 < as.Date("2021-12-01"), 1, 0)

###### Case vs non-case dist'n

A <- with(subset(ensemble.dat.immuno, Region == 0), ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
t0 <- 53

# pdf(file = "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/Immuno-Vioplots-Case-US.pdf", width = 9, height = 7)
pdf(file = paste0(spath, "/reports/Immuno-Vioplots-Case-US_",  Sys.Date(), ".pdf"))

#BA.4.5 (Severe-critical; Hybrid)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer (Hybrid; Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 1), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 1), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case& A == 1],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1], na.rm = TRUE), 1)))

#BA.4.5 (Mod-to-Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer (Hybrid; Moderate to Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 1), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 1), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case& A == 1],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1], na.rm = TRUE), 1)))

#BA.4.5 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer (Hybrid; Infection, Any Symptomatology)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 1), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 1), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case& A == 1],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 1], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 1], na.rm = TRUE), 1)))



#B1 (Severe-critical; Hybrid)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer (Hybrid; Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 1), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 1), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case& A == 1],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1], na.rm = TRUE), 1)))

#B.1 (Mod-to-Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer (Hybrid; Moderate to Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 1), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 1), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case& A == 1],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1], na.rm = TRUE), 1)))

#B.1 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer (Hybrid; Infection, Any Symptomatology)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 1), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 1), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case& A == 1],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 1], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 1),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 1], na.rm = TRUE), 1)))


##### Vaccine
#BA.4.5 (Severe-critical; Vaccine)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer (Vaccine; Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 0), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 0), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case& A == 0],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0], na.rm = TRUE), 1)))

#BA.4.5 (Mod-to-Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer (Vaccine; Moderate to Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 0), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 0), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case& A == 0],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0], na.rm = TRUE), 1)))

#BA.4.5 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 BA.4/5 Titer (Vaccine; Infection, Any Symptomatology)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 0), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 0), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case& A == 0],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[case & A == 0], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5[non.case & A == 0], na.rm = TRUE), 1)))



#B1 (Severe-critical; Vaccine)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer (Vaccine; Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 0), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 0), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case& A == 0],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0], na.rm = TRUE), 1)))

#B.1 (Mod-to-Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer (Vaccine; Moderate to Severe-critical)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 0), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 0), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case& A == 0],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0], na.rm = TRUE), 1)))

#B.1 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0
vioplot(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0],
        subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0],
        drawRect = FALSE,
        xlab = "", ylab = "nAb-ID50 Titer (AU/mL)",
        main = "BD-28 nAb-ID50 B.1 Titer (Vaccine; Infection, Any Symptomatology)", border = c("salmon", "blue"),
        col = c("white", "white"),
        lwd = 2,
        cex.axis = 1, ylim = c(0,5.4), yaxt = "n", xaxt = "n")
axis(at = c(1, 2, 3, 4), side = 2,
     labels = c(TeX("$10^1$"), TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$")))
axis(at = c(1, 2), side = 1, labels = c("Case\n(within 60 Days post BD-28)", "Non-case\n(through 60 Days post BD-28)"), tick = FALSE)

points(x = rnorm(sum(case & A == 0), mean = 1, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0],
       cex = .4, col = "black")
points(x = rnorm(sum(non.case & A == 0), mean = 2, sd = .1),
       y = subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case& A == 0],
       cex = .4, col = "black")


ppr.cutoff <- log10(2.612/2)
axis(at = log10(2.612), labels = "LOD", side = 2, las = 2)

abline(h = log10(2.612), lwd = 1.5, col = "grey", lty = 2)

text(x = 1, y = 5, paste0("N = ", sum(case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[case & A == 0], na.rm = TRUE), 1)))
text(x = 2, y = 5, paste0("N = ", sum(non.case & A == 0),
                          "\nPR = ", round(100 * mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0] >
                                                        ppr.cutoff, na.rm = TRUE),
                                           digits = 1), "%",
                          "\nGM = ", round(10^mean(subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1[non.case & A == 0], na.rm = TRUE), 1)))



dev.off()


### Inverse probability of sampling weighted means
tab.hybrid <- matrix(NA, nrow = 6, ncol = 6)
tab.hybrid[,1] <- c("BA.4/5", "BA.4/5", "BA.4/5", "B.1", "B.1", "B.1")
tab.hybrid[,2] <- c("Severe-critical", "Moderate to Severe-critical", "Infection",
                    "Severe-critical", "Moderate to Severe-critical", "Infection")

A <- with(subset(ensemble.dat.immuno, Region == 0), ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))

#### Hybrid

X <- subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5
#BA.45 (Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & non.case])

case.est <- mean(X.case[A == 1 & case])
case.ci.l <- mean(X.case[A == 1 & case]) - qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))
case.ci.u <- mean(X.case[A == 1 & case]) + qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))

non.case.est <- mean(X.non.case[A == 1 & non.case])
non.case.ci.l <- mean(X.non.case[A == 1 & non.case]) - qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))
non.case.ci.u <- mean(X.non.case[A == 1 & non.case]) + qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))

diff.est <- mean(X.non.case[A == 1 & non.case]) - mean(X.case[A == 1 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
                              var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case)),
                df = 1, lower.tail = FALSE)

tab.hybrid[1,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                              round(10^non.case.ci.l, digits = 1), ", ",
                              round(10^non.case.ci.u, digits = 1), ")"),
                       paste0(round(10^case.est, 1), " (", 
                              round(10^case.ci.l, digits = 1), ", ",
                              round(10^case.ci.u, digits = 1), ")"),
                       paste0(round(10^diff.est, 1), " (", 
                              round(10^diff.ci.l, digits = 1), ", ",
                              round(10^diff.ci.u, digits = 1), ")"),
                       signif(p.val, 6))

#BA.45 (ModtoSev)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & non.case])

case.est <- mean(X.case[A == 1 & case])
case.ci.l <- mean(X.case[A == 1 & case]) - qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))
case.ci.u <- mean(X.case[A == 1 & case]) + qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))

non.case.est <- mean(X.non.case[A == 1 & non.case])
non.case.ci.l <- mean(X.non.case[A == 1 & non.case]) - qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))
non.case.ci.u <- mean(X.non.case[A == 1 & non.case]) + qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))

diff.est <- mean(X.non.case[A == 1 & non.case]) - mean(X.case[A == 1 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
                              var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case)),
                df = 1, lower.tail = FALSE)

tab.hybrid[2,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                              round(10^non.case.ci.l, digits = 1), ", ",
                              round(10^non.case.ci.u, digits = 1), ")"),
                       paste0(round(10^case.est, 1), " (", 
                              round(10^case.ci.l, digits = 1), ", ",
                              round(10^case.ci.u, digits = 1), ")"),
                       paste0(round(10^diff.est, 1), " (", 
                              round(10^diff.ci.l, digits = 1), ", ",
                              round(10^diff.ci.u, digits = 1), ")"),
                       signif(p.val, 6))

#BA.45 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & non.case])

case.est <- mean(X.case[A == 1 & case])
case.ci.l <- mean(X.case[A == 1 & case]) - qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))
case.ci.u <- mean(X.case[A == 1 & case]) + qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))

non.case.est <- mean(X.non.case[A == 1 & non.case])
non.case.ci.l <- mean(X.non.case[A == 1 & non.case]) - qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))
non.case.ci.u <- mean(X.non.case[A == 1 & non.case]) + qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))

diff.est <- mean(X.non.case[A == 1 & non.case]) - mean(X.case[A == 1 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
                              var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case)),
                df = 1, lower.tail = FALSE)

tab.hybrid[3,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                              round(10^non.case.ci.l, digits = 1), ", ",
                              round(10^non.case.ci.u, digits = 1), ")"),
                       paste0(round(10^case.est, 1), " (", 
                              round(10^case.ci.l, digits = 1), ", ",
                              round(10^case.ci.u, digits = 1), ")"),
                       paste0(round(10^diff.est, 1), " (", 
                              round(10^diff.ci.l, digits = 1), ", ",
                              round(10^diff.ci.u, digits = 1), ")"),
                       signif(p.val, 6))



X <- subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1
#B.1 (Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & non.case])

case.est <- mean(X.case[A == 1 & case])
case.ci.l <- mean(X.case[A == 1 & case]) - qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))
case.ci.u <- mean(X.case[A == 1 & case]) + qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))

non.case.est <- mean(X.non.case[A == 1 & non.case])
non.case.ci.l <- mean(X.non.case[A == 1 & non.case]) - qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))
non.case.ci.u <- mean(X.non.case[A == 1 & non.case]) + qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))

diff.est <- mean(X.non.case[A == 1 & non.case]) - mean(X.case[A == 1 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
                              var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case)),
                df = 1, lower.tail = FALSE)

tab.hybrid[4,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                              round(10^non.case.ci.l, digits = 1), ", ",
                              round(10^non.case.ci.u, digits = 1), ")"),
                       paste0(round(10^case.est, 1), " (", 
                              round(10^case.ci.l, digits = 1), ", ",
                              round(10^case.ci.u, digits = 1), ")"),
                       paste0(round(10^diff.est, 1), " (", 
                              round(10^diff.ci.l, digits = 1), ", ",
                              round(10^diff.ci.u, digits = 1), ")"),
                       signif(p.val, 6))

#B.1 (ModtoSev)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & non.case])

case.est <- mean(X.case[A == 1 & case])
case.ci.l <- mean(X.case[A == 1 & case]) - qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))
case.ci.u <- mean(X.case[A == 1 & case]) + qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))

non.case.est <- mean(X.non.case[A == 1 & non.case])
non.case.ci.l <- mean(X.non.case[A == 1 & non.case]) - qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))
non.case.ci.u <- mean(X.non.case[A == 1 & non.case]) + qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))

diff.est <- mean(X.non.case[A == 1 & non.case]) - mean(X.case[A == 1 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
                              var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case)),
                df = 1, lower.tail = FALSE)

tab.hybrid[5,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                              round(10^non.case.ci.l, digits = 1), ", ",
                              round(10^non.case.ci.u, digits = 1), ")"),
                       paste0(round(10^case.est, 1), " (", 
                              round(10^case.ci.l, digits = 1), ", ",
                              round(10^case.ci.u, digits = 1), ")"),
                       paste0(round(10^diff.est, 1), " (", 
                              round(10^diff.ci.l, digits = 1), ", ",
                              round(10^diff.ci.u, digits = 1), ")"),
                       signif(p.val, 6))

#B.1 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 1 & non.case])

case.est <- mean(X.case[A == 1 & case])
case.ci.l <- mean(X.case[A == 1 & case]) - qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))
case.ci.u <- mean(X.case[A == 1 & case]) + qnorm(.975) * sd(X.case[A == 1 & case])/sqrt(sum(A == 1 & case))

non.case.est <- mean(X.non.case[A == 1 & non.case])
non.case.ci.l <- mean(X.non.case[A == 1 & non.case]) - qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))
non.case.ci.u <- mean(X.non.case[A == 1 & non.case]) + qnorm(.975) * sd(X.non.case[A == 1 & non.case])/sqrt(sum(A == 1 & non.case))

diff.est <- mean(X.non.case[A == 1 & non.case]) - mean(X.case[A == 1 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
         var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 1 & case])/sum(A == 1 & case) + 
                              var(X.non.case[A == 1 & non.case])/sum(A == 1 & non.case)),
                df = 1, lower.tail = FALSE)

tab.hybrid[6,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                              round(10^non.case.ci.l, digits = 1), ", ",
                              round(10^non.case.ci.u, digits = 1), ")"),
                       paste0(round(10^case.est, 1), " (", 
                              round(10^case.ci.l, digits = 1), ", ",
                              round(10^case.ci.u, digits = 1), ")"),
                       paste0(round(10^diff.est, 1), " (", 
                              round(10^diff.ci.l, digits = 1), ", ",
                              round(10^diff.ci.u, digits = 1), ")"),
                       signif(p.val, 6))


#### Vaccine
tab.vaccine <- matrix(NA, nrow = 6, ncol = 6)
tab.vaccine[,1] <- c("BA.4/5", "BA.4/5", "BA.4/5", "B.1", "B.1", "B.1")
tab.vaccine[,2] <- c("Severe-critical", "Moderate to Severe-critical", "Infection",
                     "Severe-critical", "Moderate to Severe-critical", "Infection")

A <- with(subset(ensemble.dat.immuno, Region == 0), ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))

X <- subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_BA.4.5
#BA.45 (Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & non.case])

case.est <- mean(X.case[A == 0 & case])
case.ci.l <- mean(X.case[A == 0 & case]) - qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))
case.ci.u <- mean(X.case[A == 0 & case]) + qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))

non.case.est <- mean(X.non.case[A == 0 & non.case])
non.case.ci.l <- mean(X.non.case[A == 0 & non.case]) - qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))
non.case.ci.u <- mean(X.non.case[A == 0 & non.case]) + qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))

diff.est <- mean(X.non.case[A == 0 & non.case]) - mean(X.case[A == 0 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
                              var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case)),
                df = 1, lower.tail = FALSE)

tab.vaccine[1,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                               round(10^non.case.ci.l, digits = 1), ", ",
                               round(10^non.case.ci.u, digits = 1), ")"),
                        paste0(round(10^case.est, 1), " (", 
                               round(10^case.ci.l, digits = 1), ", ",
                               round(10^case.ci.u, digits = 1), ")"),
                        paste0(round(10^diff.est, 1), " (", 
                               round(10^diff.ci.l, digits = 1), ", ",
                               round(10^diff.ci.u, digits = 1), ")"),
                        signif(p.val, 6))

#BA.45 (ModtoSev)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & non.case])

case.est <- mean(X.case[A == 0 & case])
case.ci.l <- mean(X.case[A == 0 & case]) - qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))
case.ci.u <- mean(X.case[A == 0 & case]) + qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))

non.case.est <- mean(X.non.case[A == 0 & non.case])
non.case.ci.l <- mean(X.non.case[A == 0 & non.case]) - qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))
non.case.ci.u <- mean(X.non.case[A == 0 & non.case]) + qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))

diff.est <- mean(X.non.case[A == 0 & non.case]) - mean(X.case[A == 0 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
                              var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case)),
                df = 1, lower.tail = FALSE)

tab.vaccine[2,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                               round(10^non.case.ci.l, digits = 1), ", ",
                               round(10^non.case.ci.u, digits = 1), ")"),
                        paste0(round(10^case.est, 1), " (", 
                               round(10^case.ci.l, digits = 1), ", ",
                               round(10^case.ci.u, digits = 1), ")"),
                        paste0(round(10^diff.est, 1), " (", 
                               round(10^diff.ci.l, digits = 1), ", ",
                               round(10^diff.ci.u, digits = 1), ")"),
                        signif(p.val, 6))

#BA.45 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & non.case])

case.est <- mean(X.case[A == 0 & case])
case.ci.l <- mean(X.case[A == 0 & case]) - qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))
case.ci.u <- mean(X.case[A == 0 & case]) + qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))

non.case.est <- mean(X.non.case[A == 0 & non.case])
non.case.ci.l <- mean(X.non.case[A == 0 & non.case]) - qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))
non.case.ci.u <- mean(X.non.case[A == 0 & non.case]) + qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))

diff.est <- mean(X.non.case[A == 0 & non.case]) - mean(X.case[A == 0 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
                              var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case)),
                df = 1, lower.tail = FALSE)

tab.vaccine[3,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                               round(10^non.case.ci.l, digits = 1), ", ",
                               round(10^non.case.ci.u, digits = 1), ")"),
                        paste0(round(10^case.est, 1), " (", 
                               round(10^case.ci.l, digits = 1), ", ",
                               round(10^case.ci.u, digits = 1), ")"),
                        paste0(round(10^diff.est, 1), " (", 
                               round(10^diff.ci.l, digits = 1), ", ",
                               round(10^diff.ci.u, digits = 1), ")"),
                        signif(p.val, 6))



X <- subset(ensemble.dat.immuno, Region == 0)$Year1_28dayspseudoneutid50_B.1
#B.1 (Severe-critical)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & non.case])

case.est <- mean(X.case[A == 0 & case])
case.ci.l <- mean(X.case[A == 0 & case]) - qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))
case.ci.u <- mean(X.case[A == 0 & case]) + qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))

non.case.est <- mean(X.non.case[A == 0 & non.case])
non.case.ci.l <- mean(X.non.case[A == 0 & non.case]) - qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))
non.case.ci.u <- mean(X.non.case[A == 0 & non.case]) + qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))

diff.est <- mean(X.non.case[A == 0 & non.case]) - mean(X.case[A == 0 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
                              var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case)),
                df = 1, lower.tail = FALSE)

tab.vaccine[4,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                               round(10^non.case.ci.l, digits = 1), ", ",
                               round(10^non.case.ci.u, digits = 1), ")"),
                        paste0(round(10^case.est, 1), " (", 
                               round(10^case.ci.l, digits = 1), ", ",
                               round(10^case.ci.u, digits = 1), ")"),
                        paste0(round(10^diff.est, 1), " (", 
                               round(10^diff.ci.l, digits = 1), ", ",
                               round(10^diff.ci.u, digits = 1), ")"),
                        signif(p.val, 6))

#B.1 (ModtoSev)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & non.case])

case.est <- mean(X.case[A == 0 & case])
case.ci.l <- mean(X.case[A == 0 & case]) - qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))
case.ci.u <- mean(X.case[A == 0 & case]) + qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))

non.case.est <- mean(X.non.case[A == 0 & non.case])
non.case.ci.l <- mean(X.non.case[A == 0 & non.case]) - qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))
non.case.ci.u <- mean(X.non.case[A == 0 & non.case]) + qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))

diff.est <- mean(X.non.case[A == 0 & non.case]) - mean(X.case[A == 0 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
                              var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case)),
                df = 1, lower.tail = FALSE)

tab.vaccine[5,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                               round(10^non.case.ci.l, digits = 1), ", ",
                               round(10^non.case.ci.u, digits = 1), ")"),
                        paste0(round(10^case.est, 1), " (", 
                               round(10^case.ci.l, digits = 1), ", ",
                               round(10^case.ci.u, digits = 1), ")"),
                        paste0(round(10^diff.est, 1), " (", 
                               round(10^diff.ci.l, digits = 1), ", ",
                               round(10^diff.ci.u, digits = 1), ")"),
                        signif(p.val, 6))

#B.1 (Any)
case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= t0
non.case <- subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
  subset(ensemble.dat.immuno, Region == 0)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > t0

X.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & case])
X.non.case <- X * subset(ensemble.dat.immuno, Region == 0)$wt.immuno/mean(subset(ensemble.dat.immuno, Region == 0)$wt.immuno[A == 0 & non.case])

case.est <- mean(X.case[A == 0 & case])
case.ci.l <- mean(X.case[A == 0 & case]) - qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))
case.ci.u <- mean(X.case[A == 0 & case]) + qnorm(.975) * sd(X.case[A == 0 & case])/sqrt(sum(A == 0 & case))

non.case.est <- mean(X.non.case[A == 0 & non.case])
non.case.ci.l <- mean(X.non.case[A == 0 & non.case]) - qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))
non.case.ci.u <- mean(X.non.case[A == 0 & non.case]) + qnorm(.975) * sd(X.non.case[A == 0 & non.case])/sqrt(sum(A == 0 & non.case))

diff.est <- mean(X.non.case[A == 0 & non.case]) - mean(X.case[A == 0 & case])
diff.ci.l <- diff.est - qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
diff.ci.u <- diff.est + qnorm(.975) * 
  sqrt(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
         var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case))
p.val <- pchisq(diff.est^2/(var(X.case[A == 0 & case])/sum(A == 0 & case) + 
                              var(X.non.case[A == 0 & non.case])/sum(A == 0 & non.case)),
                df = 1, lower.tail = FALSE)

tab.vaccine[6,3:6] <- c(paste0(round(10^non.case.est, 1), " (", 
                               round(10^non.case.ci.l, digits = 1), ", ",
                               round(10^non.case.ci.u, digits = 1), ")"),
                        paste0(round(10^case.est, 1), " (", 
                               round(10^case.ci.l, digits = 1), ", ",
                               round(10^case.ci.u, digits = 1), ")"),
                        paste0(round(10^diff.est, 1), " (", 
                               round(10^diff.ci.l, digits = 1), ", ",
                               round(10^diff.ci.u, digits = 1), ")"),
                        signif(p.val, 6))

tab <- cbind(c(rep("Hybrid", 6), rep("Vaccine", 6)), rbind(tab.hybrid, tab.vaccine))
colnames(tab) <- c("Group", "Antigen", "Endpoint","Non-case GM", "Case GM", "GMR", "P-value")
# write.csv(tab, "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/correlates-weighted-GM-US.csv")
write.csv(tab, paste0(spath, "/reports/correlates-weighted-GM-LatinAmerica-US_", Sys.Date(), ".csv"))

