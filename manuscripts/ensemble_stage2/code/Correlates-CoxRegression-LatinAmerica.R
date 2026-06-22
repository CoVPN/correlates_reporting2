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

# ensemble.dat <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv")
ensemble.dat <- ensemble.dat[ensemble.dat$Cohort %in% c("Cohort 1", "Cohort 2") & ensemble.dat$History %in% c("History a", "History b"),]
ensemble.dat <- ensemble.dat[ensemble.dat$BOOSTDT < as.Date("2022-05-22"),]
ensemble.dat$White <- with(ensemble.dat, ifelse(Black == 0 & Asian == 0 & NatAmer == 0 & PacIsl == 0 &
                                                  Multiracial == 0 & Notreported == 0 & Unknown == 0, 1, 0))
ensemble.dat$OtherRace <- with(ensemble.dat, ifelse(Notreported == 1 | Unknown == 1, 1, 0))

# ensemble.dat.immuno <- read.csv("/Volumes/trials/covpn/p3003/analysis/mapping_immune_correlates/post unblinded part B stage 2/adata/COVID_ENSEMBLE_stage2_mapped_20251203.csv")
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

###### Cox regression correlates analysis

##B.1
tab.hy <- matrix(NA, nrow = 2, ncol = 2)
tab.vac <- matrix(NA, nrow = 2, ncol= 2)
tab.int <- matrix(NA, nrow = 2, ncol = 2)

X <-  subset(subset(ensemble.dat.immuno, Region == 1))$Year1_28dayspseudoneutid50_B.1
A <- with(subset(ensemble.dat.immuno, Region == 1), ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
W <- with(subset(ensemble.dat.immuno, Region == 1), cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))
cal.time.num <- as.numeric(as.Date(subset(ensemble.dat.immuno, Region == 1)$SVSTDTC_YEAR.1_28.DAYS))

# complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
# 
# cox.sev.hy <- coxph(Surv(cal.time.num[A == 1 & complete], 
#                          cal.time.num[A == 1 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 1 & complete],
#                          subset(ensemble.dat.immuno, Region == 1)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 1 & complete]) ~
#                       X[A == 1 & complete] + W[A == 1 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 1 & complete],
#                     robust = TRUE, id = 1:sum(A == 1 & complete))
# cox.sev.vac <- coxph(Surv(cal.time.num[A == 0 & complete], 
#                           cal.time.num[A == 0 & complete] +  subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 0 & complete],
#                           event = subset(ensemble.dat.immuno, Region == 1)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 0 & complete]) ~
#                        X[A == 0 & complete] + W[A == 0 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 0 & complete],
#                      robust = TRUE, id = 1:sum(A == 0 & complete))
# cox.sev.int <- coxph(Surv(cal.time.num[complete], 
#                           cal.time.num[complete] +  subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[complete],
#                           event = subset(ensemble.dat.immuno, Region == 1)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[complete]) ~
#                        X[complete] * A[complete] + W[complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[complete],
#                      robust = TRUE,  id = 1:sum(complete))

complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0

cox.modtosev.hy <- coxph(Surv(cal.time.num[A == 1 & complete], 
                              cal.time.num[A == 1 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 1 & complete],
                              event = subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 1 & complete]) ~
                           X[A == 1 & complete] + W[A == 1 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 1 & complete],
                         robust = TRUE,  id = 1:sum(A == 1 & complete))
cox.modtosev.vac <- coxph(Surv(cal.time.num[A == 0 & complete], 
                               cal.time.num[A == 0 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 0 & complete],
                               event = subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 0 & complete]) ~
                            X[A == 0 & complete] + W[A == 0 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 0 & complete],
                          robust = TRUE,  id = 1:sum(A == 0 & complete))
cox.modtosev.int <- coxph(Surv(cal.time.num[complete], 
                               cal.time.num[complete] + subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[complete],
                               event = subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[complete]) ~
                            X[complete] * A[complete] + W[complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[complete],
                          robust = TRUE,  id = 1:sum(complete))

complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 0

cox.any.hy <- coxph(Surv(cal.time.num[A == 1 & complete], 
                         cal.time.num[A == 1 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[A == 1 & complete],
                         event = subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[A == 1 & complete]) ~
                      X[A == 1 & complete] + W[A == 1 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 1 & complete],
                    robust = TRUE,  id = 1:sum(A == 1 & complete))
cox.any.vac <- coxph(Surv(cal.time.num[A == 0 & complete], 
                          cal.time.num[A == 0 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[A == 0 & complete],
                          event = subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[A == 0 & complete]) ~
                       X[A == 0 & complete] + W[A == 0 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 0 & complete],
                     robust = TRUE,  id = 1:sum(A == 0 & complete))
cox.any.int <- coxph(Surv(cal.time.num[complete], 
                          cal.time.num[complete] + subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[complete],
                          event = subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[complete]) ~
                       X[complete] * A[complete] + W[complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[complete],
                     robust = TRUE,  id = 1:sum(complete))


tab.hy[1,1] <- paste0(round(summary(cox.modtosev.hy)$coefficients[1,2], 2), " (",
                      round(exp(summary(cox.modtosev.hy)$coefficients[1,1] - 
                                  qnorm(.975) *summary(cox.modtosev.hy)$coefficients[1,4]), 2),
                      ", ",
                      round(exp(summary(cox.modtosev.hy)$coefficients[1,1] +
                                  qnorm(.975) *summary(cox.modtosev.hy)$coefficients[1,4]), 2), ")")
tab.hy[1,2] <- signif(summary(cox.modtosev.hy)$coefficients[1,6])

tab.hy[2,1] <- paste0(round(summary(cox.any.hy)$coefficients[1,2], 2), " (",
                      round(exp(summary(cox.any.hy)$coefficients[1,1] - 
                                  qnorm(.975) *summary(cox.any.hy)$coefficients[1,4]), 2),
                      ", ",
                      round(exp(summary(cox.any.hy)$coefficients[1,1] +
                                  qnorm(.975) *summary(cox.any.hy)$coefficients[1,4]), 2), ")")
tab.hy[2,2] <- signif(summary(cox.any.hy)$coefficients[1,6])


tab.vac[1,1] <- paste0(round(summary(cox.modtosev.vac)$coefficients[1,2], 2), " (",
                       round(exp(summary(cox.modtosev.vac)$coefficients[1,1] - 
                                   qnorm(.975) *summary(cox.modtosev.vac)$coefficients[1,4]), 2),
                       ", ",
                       round(exp(summary(cox.modtosev.vac)$coefficients[1,1] +
                                   qnorm(.975) *summary(cox.modtosev.vac)$coefficients[1,4]), 2), ")")
tab.vac[1,2] <- signif(summary(cox.modtosev.vac)$coefficients[1,6])

tab.vac[2,1] <- paste0(round(summary(cox.any.vac)$coefficients[1,2], 2), " (",
                       round(exp(summary(cox.any.vac)$coefficients[1,1] - 
                                   qnorm(.975) *summary(cox.any.vac)$coefficients[1,4]), 2),
                       ", ",
                       round(exp(summary(cox.any.vac)$coefficients[1,1] +
                                   qnorm(.975) *summary(cox.any.vac)$coefficients[1,4]), 2), ")")
tab.vac[2,2] <- signif(summary(cox.any.vac)$coefficients[1,6])


tab.int[1,1] <- paste0(round(summary(cox.modtosev.int)$coefficients[7,2], 2), " (",
                       round(exp(summary(cox.modtosev.int)$coefficients[7,1] - 
                                   qnorm(.975) *summary(cox.modtosev.int)$coefficients[7,4]), 2),
                       ", ",
                       round(exp(summary(cox.modtosev.int)$coefficients[7,1] +
                                   qnorm(.975) *summary(cox.modtosev.int)$coefficients[7,4]), 2), ")")
tab.int[1,2] <- signif(summary(cox.modtosev.int)$coefficients[7,6])

tab.int[2,1] <- paste0(round(summary(cox.any.int)$coefficients[7,2], 2), " (",
                       round(exp(summary(cox.any.int)$coefficients[7,1] - 
                                   qnorm(.975) *summary(cox.any.int)$coefficients[7,4]), 2),
                       ", ",
                       round(exp(summary(cox.any.int)$coefficients[7,1] +
                                   qnorm(.975) *summary(cox.any.int)$coefficients[7,4]), 2), ")")
tab.int[2,2] <- signif(summary(cox.any.int)$coefficients[7,6])


tab <- cbind(tab.hy, tab.vac, tab.int)
rownames(tab) <- c("Mod to Severe-critical", "Any Infection")
colnames(tab) <- c("Hybrid (HR, CI)", "Hybrid (P-value)", 
                   "Vaccine (HR, CI)", "Vaccine (P-value)",
                   "Interaction HR (CI)", "Interaction HR (P-value)")
# write.csv(tab, "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/CoxCorrelates-B1-LatinAmerica.csv")
write.csv(tab, paste0(spath, "/reports/CoxCorrelates-B1-LatinaAmerica_", Sys.Date(), ".csv"))



##BA.4/5
tab.hy <- matrix(NA, nrow = 2, ncol = 2)
tab.vac <- matrix(NA, nrow = 2, ncol= 2)
tab.int <- matrix(NA, nrow = 2, ncol = 2)

X <-  subset(ensemble.dat.immuno, Region == 1)$Year1_28dayspseudoneutid50_BA.4.5
A <- with(subset(ensemble.dat.immuno, Region == 1), ifelse(group %in% c("Immediate Hybrid", "Delayed Hybrid"), 1, 0))
# W <- with(subset(ensemble.dat.immuno, Region == 1), cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))
W <- with(subset(ensemble.dat.immuno, Region == 1), cbind(Age, Sex, BMIDich * Sex, BMIDich * (1 - Sex)))
# W <- with(subset(ensemble.dat.immuno, Region == 1), cbind(Age, BMIDich * Sex, BMIDich * (1 - Sex),
#                                      RegionUS, RegionLA, RegionRSA, HistoryA, CalendarDateDose2Early))

# complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))})
# complete <- complete & ((subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 1 &
#                            subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc <= 90) |
#                           (subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc == 0 &
#                              subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 90))

# complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0
# 
# cox.sev.hy <- coxph(Surv(cal.time.num[A == 1 & complete], 
#                          cal.time.num[A == 1 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 1 & complete],
#                          subset(ensemble.dat.immuno, Region == 1)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 1 & complete]) ~
#                       X[A == 1 & complete] + W[A == 1 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 1 & complete],
#                     robust = TRUE, id = 1:sum(A == 1 & complete))
# cox.sev.vac <- coxph(Surv(cal.time.num[A == 0 & complete], 
#                           cal.time.num[A == 0 & complete] +  subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 0 & complete],
#                           event = subset(ensemble.dat.immuno, Region == 1)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 0 & complete]) ~
#                        X[A == 0 & complete] + W[A == 0 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 0 & complete],
#                      robust = TRUE, id = 1:sum(A == 0 & complete))
# cox.sev.int <- coxph(Surv(cal.time.num[complete], 
#                           cal.time.num[complete] +  subset(ensemble.dat.immuno, Region == 1)$AdjSevEventTimePrimaryBD35_CensorAtOutsideVacc[complete],
#                           event = subset(ensemble.dat.immuno, Region == 1)$AdjSevEventIndPrimaryBD35_CensorAtOutsideVacc[complete]) ~
#                        X[complete] * A[complete] + W[complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[complete],
#                      robust = TRUE,  id = 1:sum(complete))

complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc > 0

cox.modtosev.hy <- coxph(Surv(cal.time.num[A == 1 & complete], 
                              cal.time.num[A == 1 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 1 & complete],
                              event = subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 1 & complete]) ~
                           X[A == 1 & complete] + W[A == 1 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 1 & complete],
                         robust = TRUE,  id = 1:sum(A == 1 & complete))
cox.modtosev.vac <- coxph(Surv(cal.time.num[A == 0 & complete], 
                               cal.time.num[A == 0 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[A == 0 & complete],
                               event = subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[A == 0 & complete]) ~
                            X[A == 0 & complete] + W[A == 0 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 0 & complete],
                          robust = TRUE,  id = 1:sum(A == 0 & complete))
cox.modtosev.int <- coxph(Surv(cal.time.num[complete], 
                               cal.time.num[complete] + subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventTimePrimaryBD35_CensorAtOutsideVacc[complete],
                               event = subset(ensemble.dat.immuno, Region == 1)$AdjModtoSevEventIndPrimaryBD35_CensorAtOutsideVacc[complete]) ~
                            X[complete] * A[complete] + W[complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[complete],
                          robust = TRUE,  id = 1:sum(complete))

complete <- apply(cbind(A, W, X), 1, function(x){all(!is.na(x))}) & subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc > 0

cox.any.hy <- coxph(Surv(cal.time.num[A == 1 & complete], 
                         cal.time.num[A == 1 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[A == 1 & complete],
                         event = subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[A == 1 & complete]) ~
                      X[A == 1 & complete] + W[A == 1 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 1 & complete],
                    robust = TRUE,  id = 1:sum(A == 1 & complete))
cox.any.vac <- coxph(Surv(cal.time.num[A == 0 & complete], 
                          cal.time.num[A == 0 & complete] + subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[A == 0 & complete],
                          event = subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[A == 0 & complete]) ~
                       X[A == 0 & complete] + W[A == 0 & complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[A == 0 & complete],
                     robust = TRUE,  id = 1:sum(A == 0 & complete))
cox.any.int <- coxph(Surv(cal.time.num[complete], 
                          cal.time.num[complete] + subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventTimePrimaryBD35_CensorAtOutsideVacc[complete],
                          event = subset(ensemble.dat.immuno, Region == 1)$AdjAnyEventIndPrimaryBD35_CensorAtOutsideVacc[complete]) ~
                       X[complete] * A[complete] + W[complete,], weights = subset(ensemble.dat.immuno, Region == 1)$wt.immuno[complete],
                     robust = TRUE,  id = 1:sum(complete))


tab.hy[1,1] <- paste0(round(summary(cox.modtosev.hy)$coefficients[1,2], 2), " (",
                      round(exp(summary(cox.modtosev.hy)$coefficients[1,1] - 
                                  qnorm(.975) *summary(cox.modtosev.hy)$coefficients[1,4]), 2),
                      ", ",
                      round(exp(summary(cox.modtosev.hy)$coefficients[1,1] +
                                  qnorm(.975) *summary(cox.modtosev.hy)$coefficients[1,4]), 2), ")")
tab.hy[1,2] <- signif(summary(cox.modtosev.hy)$coefficients[1,6])

tab.hy[2,1] <- paste0(round(summary(cox.any.hy)$coefficients[1,2], 2), " (",
                      round(exp(summary(cox.any.hy)$coefficients[1,1] - 
                                  qnorm(.975) *summary(cox.any.hy)$coefficients[1,4]), 2),
                      ", ",
                      round(exp(summary(cox.any.hy)$coefficients[1,1] +
                                  qnorm(.975) *summary(cox.any.hy)$coefficients[1,4]), 2), ")")
tab.hy[2,2] <- signif(summary(cox.any.hy)$coefficients[1,6])


tab.vac[1,1] <- paste0(round(summary(cox.modtosev.vac)$coefficients[1,2], 2), " (",
                       round(exp(summary(cox.modtosev.vac)$coefficients[1,1] - 
                                   qnorm(.975) *summary(cox.modtosev.vac)$coefficients[1,4]), 2),
                       ", ",
                       round(exp(summary(cox.modtosev.vac)$coefficients[1,1] +
                                   qnorm(.975) *summary(cox.modtosev.vac)$coefficients[1,4]), 2), ")")
tab.vac[1,2] <- signif(summary(cox.modtosev.vac)$coefficients[1,6])

tab.vac[2,1] <- paste0(round(summary(cox.any.vac)$coefficients[1,2], 2), " (",
                       round(exp(summary(cox.any.vac)$coefficients[1,1] - 
                                   qnorm(.975) *summary(cox.any.vac)$coefficients[1,4]), 2),
                       ", ",
                       round(exp(summary(cox.any.vac)$coefficients[1,1] +
                                   qnorm(.975) *summary(cox.any.vac)$coefficients[1,4]), 2), ")")
tab.vac[2,2] <- signif(summary(cox.any.vac)$coefficients[1,6])


tab.int[1,1] <- paste0(round(summary(cox.modtosev.int)$coefficients[7,2], 2), " (",
                       round(exp(summary(cox.modtosev.int)$coefficients[7,1] - 
                                   qnorm(.975) *summary(cox.modtosev.int)$coefficients[7,4]), 2),
                       ", ",
                       round(exp(summary(cox.modtosev.int)$coefficients[7,1] +
                                   qnorm(.975) *summary(cox.modtosev.int)$coefficients[7,4]), 2), ")")
tab.int[1,2] <- signif(summary(cox.modtosev.int)$coefficients[7,6])

tab.int[2,1] <- paste0(round(summary(cox.any.int)$coefficients[7,2], 2), " (",
                       round(exp(summary(cox.any.int)$coefficients[7,1] - 
                                   qnorm(.975) *summary(cox.any.int)$coefficients[7,4]), 2),
                       ", ",
                       round(exp(summary(cox.any.int)$coefficients[7,1] +
                                   qnorm(.975) *summary(cox.any.int)$coefficients[7,4]), 2), ")")
tab.int[2,2] <- signif(summary(cox.any.int)$coefficients[7,6])


tab <- cbind(tab.hy, tab.vac, tab.int)
rownames(tab) <- c("Mod to Severe-critical", "Any Infection")
colnames(tab) <- c("Hybrid (HR, CI)", "Hybrid (P-value)", 
                   "Vaccine (HR, CI)", "Vaccine (P-value)",
                   "Interaction HR (CI)", "Interaction HR (P-value)")
# write.csv(tab, "/Volumes/ahudson/Ensemble-stage2-correlates/EfficacyImmunoCorrelates/CoxCorrelates-BA45-LatinAmerica.csv")
write.csv(tab, paste0(spath, "/reports/CoxCorrelates-BA45-LatinAmerica_", Sys.Date(), ".csv"))
