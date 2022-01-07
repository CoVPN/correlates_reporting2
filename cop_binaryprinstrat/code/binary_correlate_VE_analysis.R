#Sys.setenv(TRIAL = "moderna_real"); Args=c(COR="D57", 1); Sys.setenv(VERBOSE = 1) 
# TRIAL: moderna_mock  moderna_real  janssen_pooled_mock  janssen_pooled_real  janssen_na_mock  hvtn705
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
    
source(here::here("..", "_common.R"))
#-----------------------------------------------


if(!exists("Args")) Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==2) {
    if(Args[2]=="1") approach2=FALSE else if (Args[2]=="2") approach2=TRUE else stop("The second arg, which is the indicator for the approach, has to be either 1 or 2.")
} else {
    stop("Two parameters are expected. The first one is COR, and the second one is 1 or 2, indicator for the approach.")
}



# path for figures and tables etc
save.results.to = here::here("output"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/"); if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))


############################################################################################################
# Principal stratification inference
# Method of Gilbert, Blette, Shepherd, Hudgens (2020, J Causal Inference)

# Bryan Blette's R package psbinary
library(psbinary)
library(xtable)
#help(analyze_NEE)

data <- dat.mock

# The "NEE" and "NEH" methods are included in the code, the user specifies one or the other.
# NEH is always the relevant one for COVID-19 vaccines.
method <- "NEH"


#############################################################
# Conduct tpeak time point marker correlates analyses 

# Define variables (in the notation of the general immune correlates suite) 
# by mapping to the Moderna correlates data set

# Note: Cannot start with subsetting on ph1.tpeak==1 ptids, because these already 
# exclude participants with an early infection before time tau (tpeak visit).  
# Need to set up the Y_tau variable for the Blette code to work properly.

Ph2ptids.tpeak <-   ifelse(data$ph2,1,0)

Trt <- data$Trt
Delta.Day1 <- data$EventIndPrimaryD1
Delta.tpeak <- data$EventIndPrimary
Ttilde.tpeak <- data$EventTimePrimary

# Define Y_tautpeak: takes value 1 if a SARS-CoV-2 infection or COVID-19 case by < tpeaklag days post tpeak visit
#Earlyinfectiontpeak <- data[["EarlyinfectionD"%.%tpeak]]
#Earlyendpointtpeak  <- data[["EarlyendpointD"%.%tpeak]]
#Y_tautpeak <- ifelse((Delta.tpeak==1 & Ttilde.tpeak < tpeaklag) | Earlyinfectiontpeak==1 | Earlyendpointtpeak==1,1,0)
# YF: it may be better to define this variable as follows
Y_tautpeak <- 1-data$ph1

W1 <- data$wt
W1[Ph2ptids.tpeak==0] <- 0

#################################################################################
# Analyze the set of biomarkers in the vector markersforanalysis (names in the 'data' data frame)

markersforanalysis <- paste0("Day", tpeak, assays)

# UPDATE
# The program can do the data analysis in two ways.  
#
# Approach 1: With binary markers defined by cutting each quantitative
# marker at each of a fixed set of percentiles analysispercs. The same analysispercs
# are applied across all of the antibody markers (e.g., IgG Spike, PsV nAb ID50, etc.) 
#
# Approach 2: With binary markers defined by a single cut-point for each antibody marker.  The former 
# is useful for studying results varying over percentile cut-points, and the latter for focusing on a 
# special cut-point of interest such as below vs. above the positivity cut-off (bAb) or below vs. above 
# the LOD (neutralization and ADCP markers).
#
# The Boolean approach2 defines which approach is used.

# Define the lower limits of the antibody markers (in IU), which are used for both Approach 1 and Approach 2
# at least for plotting.
#lowerlimitsassays <- c(0.3076,1.593648,0.242,1.502)
# Youyi: I need help with the line above, to make it grab from common.R the appropriate lower cutoffs for the 4 assays 
#        used in the Moderna correlates (IgG Spike, IgG RBD, PsV nAb ID50, PsV nAb ID80).  And then similarly to grab the 
#        appropriate lower cutoffs for the 4 assays used in the ENSEMBLE correlates (IgG Spike, IgG RBD, PsV nAb ID50, 
#        ADCP).   
lowerlimitsassays <- llods[assays]

#
# Currently the code is set up to use Approach 1.  Changing approach2 to TRUE would implement the code using
# Approach 2.


# Approach 1:
if(!approach2) {
analysispercs <- c(0.4,0.5,0.6) }

# Approach 2:
if(approach2) {
analysispercs <- 0.5 # Under Approach 2, the only info needed is that length(analysispercs)==1
approach2cutpoints <- log10(lowerlimitsassays)
}
# END UPDATE

# Conduct analyses for all markers in markersforanalysis for each dichotomization in analysispercs
# The following parameters are used in plotting as well as in the inferential analyses

# UPDATE
labelsmarkersforanalysis <- paste0("D",tpeak,assays,".",rep(round(10*analysispercs),each=length(assays)))
# END UPDATE

#labelsmarkersforanalysis <- c("D57Spike.4","D57RBD.4","D57ID50.4","D57ID80.4",
#                              "D57Spike.5","D57RBD.5","D57ID50.5","D57ID80.5",
#                              "D57Spike.6","D57RBD.6","D57ID50.6","D57ID80.6")

# UPDATE
if(approach2) { labelsmarkersforanalysis <- paste0("D",tpeak,assays,".",rep(0,each=length(assays))) }
# END UPDATE

#axislabelsmarkers <- c("Day 57 IgG Spike (IU/ml)","Day 57 IgG RBD (IU/ml)","Day 57 nAb titer (IU50/ml)",
#                       "Day 57 nAb titer (IU80/ml)")
axislabelsmarkers <- paste0("Day ", tpeak, " ", assay_labels_short)

################################################################################
# Generic code for calling Bryan Blette's R package

# UPDATE
# MM is the number of markers considered (e.g., IgG Spike, IgG RBD, PsV nAb ID50, PsV nAb ID80)
# END UPDATE
MM <- length(markersforanalysis)

datareportmarkerdata <- vector("list", MM*length(analysispercs))

cutpoints <- matrix(rep(NA,MM*length(analysispercs)),nrow=MM)

for (ii in 1:MM) {
for (jj in 1:length(analysispercs)) {
    
    R <- ifelse(!is.na(data[,markersforanalysis[ii]]),1,0)
    R[is.na(R)] <- 0
    R[Ph2ptids.tpeak==0] <- 0
    
    S_star <- data[,markersforanalysis[ii]]
    
    #UPDATE
    if(!approach2) {
    cutpt <- quantile(S_star[!is.na(S_star) & Delta.Day1==0],prob=analysispercs[jj])
    S_star <- ifelse(S_star > cutpt,1,0)
    S_star[is.na(S_star)] <- 0
    cutpoints[ii,jj] <- round(cutpt,4)
    }
    
    if(approach2) {
    cutpoints[ii,jj] <- approach2cutpoints[ii] }
    
    Z <- Trt
    Y <- Delta.tpeak
    Y_tau <- Y_tautpeak
    df <- data.frame(Z, S_star, R, Y_tau, Y)
    
    # Check that have correct case counts:
    cat(paste("Marker ",ii),"\n")
    cat(paste("Percentile defining cut point = ",analysispercs[jj]),"\n")
    cat(paste("Marker level cut point ",cutpoints[ii,jj]),"\n")
    print(table(S_star[R==1 & Z==1 & Y_tau==0],Y[R==1 & Z==1 & Y_tau==0]))
    
    # analyze_NEE() has an option ‘brange’ for the single sensitivity parameter
    # analyze_NEB() has ‘brange0’ and ‘brange1’ for the two sensitivity parameters
    # analyze_NEH() has ‘brange0’ up to ‘brange3’ to cover the 4 sensitivity parameters
    
    # Results on differences across the two marker subgroups are reported on the scale 
    # (1-VElow)/(1-VEhigh).
    # This is done by applying the method using constrast = "logRR" and transforming results
    # to the (1-VElow)/(1-VEhigh) scale.
    
    # Analysis under the No Early Effect (NEE) assumption:
    # Three scenarios of the degree of selection bias:
    #   None: selection bias odds ratio 1 to 1
    #        Med: selection bias odds ratio 0.75 to 1.33
    #       High: selection bias odds ratio 0.5 to 2.0
    
    if (method=="NEE") {
        
        output_S1_none <- analyze_NEE(df, brange = c(0, 0), design = "other", contrast = "logRR", weights = W1)
        output_S1_med <- analyze_NEE(df, contrast = "logRR", weights = W1,
                                       design = "other", brange = c(log(0.75),-log(0.75)))
        output_S1_high <- analyze_NEE(df, contrast = "logRR", weights = W1,
                                        design = "other", brange = c(log(0.5),-log(0.5)))
        # Create a table of NEE results:
        igintVEhigh_none <- c(1-exp(output_S1_none$CEP_10_II[2]),1-exp(output_S1_none$CEP_10_II[1]))
        igintVElow_none <- c(1-exp(output_S1_none$CEP_00_II[2]),1-exp(output_S1_none$CEP_00_II[1]))
        euiintVEhigh_none <- c(1-exp(output_S1_none$CEP_10_EUI[2]),1-exp(output_S1_none$CEP_10_EUI[1]))
        euiintVElow_none <- c(1-exp(output_S1_none$CEP_00_EUI[2]),1-exp(output_S1_none$CEP_00_EUI[1]))
        # (1 - VElow)/(1 - VEhigh) "How much greater is VE for high marker than low marker"
        igratioRRShighoverlow_none <- 1/c(exp(output_S1_none$CEP_diff_II[2]),exp(output_S1_none$CEP_diff_II[1]))
        euiratioRRShighoverlow_none <- 1/c(exp(output_S1_none$CEP_diff_EUI[2]),exp(output_S1_none$CEP_diff_EUI[1]))
        
        igintVEhigh_med <- c(1-exp(output_S1_med$CEP_10_II[2]),1-exp(output_S1_med$CEP_10_II[1]))
        igintVElow_med <- c(1-exp(output_S1_med$CEP_00_II[2]),1-exp(output_S1_med$CEP_00_II[1]))
        euiintVEhigh_med <- c(1-exp(output_S1_med$CEP_10_EUI[2]),1-exp(output_S1_med$CEP_10_EUI[1]))
        euiintVElow_med <- c(1-exp(output_S1_med$CEP_00_EUI[2]),1-exp(output_S1_med$CEP_00_EUI[1]))
        # (1 - VElow)/(1 - VEhigh) "How much greater is VE for high marker than low marker"
        igratioRRShighoverlow_med <- 1/c(exp(output_S1_med$CEP_diff_II[2]),exp(output_S1_med$CEP_diff_II[1]))
        euiratioRRShighoverlow_med <- 1/c(exp(output_S1_med$CEP_diff_EUI[2]),exp(output_S1_med$CEP_diff_EUI[1]))
        
        igintVEhigh_high <- c(1-exp(output_S1_high$CEP_10_II[2]),1-exp(output_S1_high$CEP_10_II[1]))
        igintVElow_high <- c(1-exp(output_S1_high$CEP_00_II[2]),1-exp(output_S1_high$CEP_00_II[1]))
        euiintVEhigh_high <- c(1-exp(output_S1_high$CEP_10_EUI[2]),1-exp(output_S1_high$CEP_10_EUI[1]))
        euiintVElow_high <- c(1-exp(output_S1_high$CEP_00_EUI[2]),1-exp(output_S1_high$CEP_00_EUI[1]))
        # (1 - VElow)/(1 - VEhigh) "How much greater is VE for high marker than low marker"
        igratioRRShighoverlow_high <- 1/c(exp(output_S1_high$CEP_diff_II[2]),exp(output_S1_high$CEP_diff_II[1]))
        euiratioRRShighoverlow_high <- 1/c(exp(output_S1_high$CEP_diff_EUI[2]),exp(output_S1_high$CEP_diff_EUI[1]))
        
        datareportmarkerdata[[MM*(jj-1)+ii]]   <- data.frame(marker=c(labelsmarkersforanalysis[MM*(jj-1)+ii],labelsmarkersforanalysis[MM*(jj-1)+ii],labelsmarkersforanalysis[MM*(jj-1)+ii]),
                                  Sens=c("None","Med","High"),
                                  IILowl=c(igintVElow_none[1],igintVElow_med[1],igintVElow_high[1]),
                                  IILowu=c(igintVElow_none[2],igintVElow_med[2],igintVElow_high[2]),
                                  EUILowl=c(euiintVElow_none[1],euiintVElow_med[1],euiintVElow_high[1]),
                                  EUILowu=c(euiintVElow_none[2],euiintVElow_med[2],euiintVElow_high[2]),
                                  IIHighl=c(igintVEhigh_none[1],igintVEhigh_med[1],igintVEhigh_high[1]),
                                  IIHighu=c(igintVEhigh_none[2],igintVEhigh_med[2],igintVEhigh_high[2]),
                                  EUIHighl=c(euiintVEhigh_none[1],euiintVEhigh_med[1],euiintVEhigh_high[1]),
                      EUIHighu=c(euiintVEhigh_none[2],euiintVEhigh_med[2],euiintVEhigh_high[2]),
                                  IIcontrastl=c(igratioRRShighoverlow_none[1],igratioRRShighoverlow_med[1],igratioRRShighoverlow_high[1]),
                      IIcontrastu=c(igratioRRShighoverlow_none[2],igratioRRShighoverlow_med[2],igratioRRShighoverlow_high[2]),
                                  EUIcontrastl=c(euiratioRRShighoverlow_none[1],euiratioRRShighoverlow_med[1],euiratioRRShighoverlow_high[1]),
                      EUIcontrastu=c(euiratioRRShighoverlow_none[2],euiratioRRShighoverlow_med[2],euiratioRRShighoverlow_high[2]))
        
        colnames(datareportmarkerdata[[M*(jj-1)+ii]]) <- c("Marker","Sens","Ilol","Ilou","Elol","Elou",
                                            "Ihil","Ihiu","Ehil","Ehiu","Icnl","Icnu",
                                            "Ecnl","Ecnu") 
    }
    
    
    if (method=="NEH") {
        # Prefer the No Early Harm (NEH) design given the evidence for early VE > 0%.
        
        output_S1_none <- analyze_NEH(df, contrast = "logRR", weights = W1,
                                        design = "other", brange0 = c(0, 0), brange1 = c(0,0), 
                                        brange2 = c(0,0), brange3 = c(0,0))
        output_S1_med <- analyze_NEH(df, contrast = "logRR", weights = W1,
                                       design = "other", brange0 = c(log(0.75),-log(0.75)), brange1 = c(log(0.75),-log(0.75)), 
                                        brange2 = c(log(0.75),-log(0.75)), brange3 = c(log(0.75),-log(0.75)))
        output_S1_high <- analyze_NEH(df, contrast = "logRR", weights = W1,
                                        design = "other", brange0 = c(log(0.5),-log(0.5)), brange1 = c(log(0.5),-log(0.5)), 
                                        brange2 = c(log(0.5),-log(0.5)), brange3 = c(log(0.5),-log(0.5)))
        
        # Create a table of NEH results:
        igintVEhigh_none <- c(1-exp(output_S1_none$CEP_10_II[2]),1-exp(output_S1_none$CEP_10_II[1]))
        igintVElow_none <- c(1-exp(output_S1_none$CEP_00_II[2]),1-exp(output_S1_none$CEP_00_II[1]))
        euiintVEhigh_none <- c(1-exp(output_S1_none$CEP_10_EUI[2]),1-exp(output_S1_none$CEP_10_EUI[1]))
        euiintVElow_none <- c(1-exp(output_S1_none$CEP_00_EUI[2]),1-exp(output_S1_none$CEP_00_EUI[1]))
        # (1 - VElow)/(1 - VEhigh) "How much greater is VE for high marker than low marker"
        igratioRRShighoverlow_none <- 1/c(exp(output_S1_none$CEP_diff_II[2]),exp(output_S1_none$CEP_diff_II[1]))
        euiratioRRShighoverlow_none <- 1/c(exp(output_S1_none$CEP_diff_EUI[2]),exp(output_S1_none$CEP_diff_EUI[1]))
        
        igintVEhigh_med <- c(1-exp(output_S1_med$CEP_10_II[2]),1-exp(output_S1_med$CEP_10_II[1]))
        igintVElow_med <- c(1-exp(output_S1_med$CEP_00_II[2]),1-exp(output_S1_med$CEP_00_II[1]))
        euiintVEhigh_med <- c(1-exp(output_S1_med$CEP_10_EUI[2]),1-exp(output_S1_med$CEP_10_EUI[1]))
        euiintVElow_med <- c(1-exp(output_S1_med$CEP_00_EUI[2]),1-exp(output_S1_med$CEP_00_EUI[1]))
        # (1 - VElow)/(1 - VEhigh) "How much greater is VE for high marker than low marker"
        igratioRRShighoverlow_med <- 1/c(exp(output_S1_med$CEP_diff_II[2]),exp(output_S1_med$CEP_diff_II[1]))
        euiratioRRShighoverlow_med <- 1/c(exp(output_S1_med$CEP_diff_EUI[2]),exp(output_S1_med$CEP_diff_EUI[1]))
        
        igintVEhigh_high <- c(1-exp(output_S1_high$CEP_10_II[2]),1-exp(output_S1_high$CEP_10_II[1]))
        igintVElow_high <- c(1-exp(output_S1_high$CEP_00_II[2]),1-exp(output_S1_high$CEP_00_II[1]))
        euiintVEhigh_high <- c(1-exp(output_S1_high$CEP_10_EUI[2]),1-exp(output_S1_high$CEP_10_EUI[1]))
        euiintVElow_high <- c(1-exp(output_S1_high$CEP_00_EUI[2]),1-exp(output_S1_high$CEP_00_EUI[1]))
        # (1 - VElow)/(1 - VEhigh) "How much greater is VE for high marker than low marker"
        igratioRRShighoverlow_high <- 1/c(exp(output_S1_high$CEP_diff_II[2]),exp(output_S1_high$CEP_diff_II[1]))
        euiratioRRShighoverlow_high <- 1/c(exp(output_S1_high$CEP_diff_EUI[2]),exp(output_S1_high$CEP_diff_EUI[1]))
        
        datareportmarkerdata[[MM*(jj-1)+ii]]  <- data.frame(marker=c(labelsmarkersforanalysis[MM*(jj-1)+ii],labelsmarkersforanalysis[MM*(jj-1)+ii],labelsmarkersforanalysis[MM*(jj-1)+ii]),
                                  Sens=c("None","Med","High"),
                                  IILowl=c(igintVElow_none[1],igintVElow_med[1],igintVElow_high[1]),
                                  IILowu=c(igintVElow_none[2],igintVElow_med[2],igintVElow_high[2]),
                                  EUILowl=c(euiintVElow_none[1],euiintVElow_med[1],euiintVElow_high[1]),
                                  EUILowu=c(euiintVElow_none[2],euiintVElow_med[2],euiintVElow_high[2]),
                                  IIHighl=c(igintVEhigh_none[1],igintVEhigh_med[1],igintVEhigh_high[1]),
                                  IIHighu=c(igintVEhigh_none[2],igintVEhigh_med[2],igintVEhigh_high[2]),
                                  EUIHighl=c(euiintVEhigh_none[1],euiintVEhigh_med[1],euiintVEhigh_high[1]),
                                  EUIHighu=c(euiintVEhigh_none[2],euiintVEhigh_med[2],euiintVEhigh_high[2]),
                                  IIcontrastl=c(igratioRRShighoverlow_none[1],igratioRRShighoverlow_med[1],igratioRRShighoverlow_high[1]),
                      IIcontrastu=c(igratioRRShighoverlow_none[2],igratioRRShighoverlow_med[2],igratioRRShighoverlow_high[2]),
                                  EUIcontrastl=c(euiratioRRShighoverlow_none[1],euiratioRRShighoverlow_med[1],euiratioRRShighoverlow_high[1]),
                      EUIcontrastu=c(euiratioRRShighoverlow_none[2],euiratioRRShighoverlow_med[2],euiratioRRShighoverlow_high[2]))
        
        colnames(datareportmarkerdata[[MM*(jj-1)+ii]]) <- c("Marker","Sens","Ilol","Ilou","Elol","Elou",
                                            "Ihil","Ihiu","Ehil","Ehiu","Icnl","Icnu",
                                            "Ecnl","Ecnu") 
    }
}
}

datareport <- do.call(rbind, datareportmarkerdata)

# Write out the data results to a file
write.table(datareport, file=paste0(save.results.to, "PSbinaryNEH",study_name,"resultsDay",tpeak,".csv"), row.names = F, sep=",")
#write.table(datareport, file="H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/PSbinaryNEHresultsDay57.csv", row.names = F, sep=",")


# Print out a latex table for drawing into an RMD/latex report:

comment <- list(pos = list(0), command = NULL)
comment$pos[[1]] <- c(nrow(datareport))
comment$command <- paste0("\\hline \\hline 
                          \\multicolumn{", NCOL(datareport), "}{l}{\\footnotesize * None: beta sensitivity parameters log(1.0) and -log(1.0)} \\\\\n
                          \\multicolumn{", NCOL(datareport), "}{l}{\\footnotesize \\quad Med: beta sensitivity parameters log(0.75) and -log(0.75)} \\\\\n
                          \\multicolumn{", NCOL(datareport), "}{l}{\\footnotesize \\quad High: beta sensitivity parameters log(0.5) and -log(0.5)}")

print(xtable(datareport,
caption=paste0(study_name, ": Correlates of Vaccine Efficacy Results by Gilbert et al. (2020) 
Method for High vs. Low Marker Subgroups Under No Early Harm Assumption with 
Sensitivity Analysis Scenarios*")), 
caption.placement="top",
add.to.row = comment,
hline.after = c(-1,0),
label= paste0("tab:PrincstratbinaryD",tpeak,"markers"),
file=paste0(save.results.to, "PSbinaryNEHresultsDay",tpeak,".tex"))
#file="H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/PSbinaryNEHresultsDay57.tex")

# UPDATE  
if(!approach2) {
cutpoints <- data.frame(round(10^cutpoints))
#colnames(cutpoints) <- c(paste("Perc.",analysispercs[1]),
#                         paste("Perc.",analysispercs[2]),
#                         paste("Perc.",analysispercs[3]))
# The code now accommodates an arbitrary length of analysispercs
colnamescutpoints <- rep(NA,length(analysispercs))
for (k in 1:length(analysispercs)) {
colnamescutpoints[k] <- paste("Perc.",analysispercs[k]) }
colnames(cutpoints) <- colnamescutpoints
}
if(approach2) {
#cutpoints only has one column
cutpoints <- data.frame(round(10^cutpoints))
colnames(cutpoints) <- "LLimit"
}
# END UPDATE


rownames(cutpoints) <- paste0("D", tpeak, assays)
print(xtable(cutpoints,
caption=paste0(study_name, ": Cut-points Defining High and Low Marker Subgroups")), 
caption.placement="top",
label= paste0("tab:PrincstratbinaryD",tpeak,"markerscutpoints"),
file=paste0(save.results.to, "PSbinaryNEHresultsDay",tpeak,"cutpoints.tex"))
#file="H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/PSbinaryNEHresultsDay57cutpoints.tex")


######################################################################
# Plot the results, one plot for each of the tpeak time point antibody markers

# Helper plotting function
# Plotting code for a given immune marker and cutpoint
plotprincstratbinary <- function(datareport,markernoncases,markername,xlims,xlab1,filenameplot,main1,
                        atxaxis,labelsxaxis,cutpt,
                        leftptii,leftptci,leftpteui,rightptii,rightptci,rightpteui) {
    
    IlolNone <- datareport$Ilol[datareport$Marker==markername & datareport$Sens=="None"]
    IlouNone <- datareport$Ilou[datareport$Marker==markername & datareport$Sens=="None"]
    IlolMed <- datareport$Ilol[datareport$Marker==markername & datareport$Sens=="Med"]
    IlouMed <- datareport$Ilou[datareport$Marker==markername & datareport$Sens=="Med"]
    IlolHigh <- datareport$Ilol[datareport$Marker==markername & datareport$Sens=="High"]
    IlouHigh <- datareport$Ilou[datareport$Marker==markername & datareport$Sens=="High"]
    
    IhilNone <- datareport$Ihil[datareport$Marker==markername & datareport$Sens=="None"]
    IhiuNone <- datareport$Ihiu[datareport$Marker==markername & datareport$Sens=="None"]
    IhilMed <- datareport$Ihil[datareport$Marker==markername & datareport$Sens=="Med"]
    IhiuMed <- datareport$Ihiu[datareport$Marker==markername & datareport$Sens=="Med"]
    IhilHigh <- datareport$Ihil[datareport$Marker==markername & datareport$Sens=="High"]
    IhiuHigh <- datareport$Ihiu[datareport$Marker==markername & datareport$Sens=="High"]
    
    ElolNone <- datareport$Elol[datareport$Marker==markername & datareport$Sens=="None"]
    ElouNone <- datareport$Elou[datareport$Marker==markername & datareport$Sens=="None"]
    ElolMed <- datareport$Elol[datareport$Marker==markername & datareport$Sens=="Med"]
    ElouMed <- datareport$Elou[datareport$Marker==markername & datareport$Sens=="Med"]
    ElolHigh <- datareport$Elol[datareport$Marker==markername & datareport$Sens=="High"]
    ElouHigh <- datareport$Elou[datareport$Marker==markername & datareport$Sens=="High"]
    
    EhilNone <- datareport$Ehil[datareport$Marker==markername & datareport$Sens=="None"]
    EhiuNone <- datareport$Ehiu[datareport$Marker==markername & datareport$Sens=="None"]
    EhilMed <- datareport$Ehil[datareport$Marker==markername & datareport$Sens=="Med"]
    EhiuMed <- datareport$Ehiu[datareport$Marker==markername & datareport$Sens=="Med"]
    EhilHigh <- datareport$Ehil[datareport$Marker==markername & datareport$Sens=="High"]
    EhiuHigh <- datareport$Ehiu[datareport$Marker==markername & datareport$Sens=="High"]
    
    RatioIlNone <- datareport$Icnl[datareport$Marker==markername & datareport$Sens=="None"]
    RatioIuNone <- datareport$Icnu[datareport$Marker==markername & datareport$Sens=="None"]
    RatioIlMed <-  datareport$Icnl[datareport$Marker==markername & datareport$Sens=="Med"]
    RatioIuMed <-  datareport$Icnu[datareport$Marker==markername & datareport$Sens=="Med"]
    RatioIlHigh <- datareport$Icnl[datareport$Marker==markername & datareport$Sens=="High"]
    RatioIuHigh <- datareport$Icnu[datareport$Marker==markername & datareport$Sens=="High"]
    
    RatioElNone <- datareport$Ecnl[datareport$Marker==markername & datareport$Sens=="None"]
    RatioEuNone <- datareport$Ecnu[datareport$Marker==markername & datareport$Sens=="None"]
    RatioElMed <-  datareport$Ecnl[datareport$Marker==markername & datareport$Sens=="Med"]
    RatioEuMed <-  datareport$Ecnu[datareport$Marker==markername & datareport$Sens=="Med"]
    RatioElHigh <- datareport$Ecnl[datareport$Marker==markername & datareport$Sens=="High"]
    RatioEuHigh <- datareport$Ecnu[datareport$Marker==markername & datareport$Sens=="High"]
    
    d <- density(markernoncases,na.rm=T)
    #fromv <- min(data$Day29pseudoneutid50[data$EventIndPrimaryD1==0],na.rm=T)
    #tov <- max(data$Day29pseudoneutid50[data$EventIndPrimaryD1==0],na.rm=T)
    
    pdf(filenameplot) 
    par(cex.axis=1.3,cex.lab=1.3,cex.main=1.3,oma=c(6,1,0,0),las=1)
    plot(d,xlim=xlims,ylim=c(0,1),axes=FALSE,
    xlab=xlab1,ylab="Vaccine Efficacy (%)",main=main1) 
    polygon(d, col="gold", border="purple")
    axis(1,at=atxaxis,labels=labelsxaxis) 
    axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,20,40,60,80,100))
    abline(v=cutpt,lty=5,lwd=2,col="purple")
    
    #leftpt <- (min(data$Day29pseudoneutid50[data$EventIndPrimaryD1==0],na.rm=T)+log10(cutpoints[3,3]))/2
    #rightpt <- (max(data$Day29pseudoneutid50[data$EventIndPrimaryD1==0],na.rm=T)+log10(cutpoints[3,3]))/2
    
    points(leftptii,IlolNone)
    arrows(leftptii,IlolMed,leftptii,IlouMed,length=0.05,angle=90,code=3,lwd=2)
    #arrows(leftptii,IlolHigh,leftptii,IlouHigh,length=0.05,angle=90,code=3,lty=2,lwd=2)
    text(leftptii+0.3,IlolNone,round(100*IlolNone))
    text(leftptii+0.3,IlolMed-0.02,round(100*IlolMed))
    text(leftptii+0.3,IlouMed+0.02,round(100*IlouMed))
    #text(leftptii+0.3,IlolHigh-0.02,round(100*IlolHigh))
    #text(leftptii+0.3,IlouHigh+0.02,round(100*IlouHigh))
    
    arrows(leftptci,ElolNone,leftptci,ElouNone,length=0.05,angle=90,code=3,lwd=2,lty=2,col="brown")
    text(leftptci+0.3,ElolNone-0.02,round(100*ElolNone))
    text(leftptci+0.3,ElouNone+0.02,round(100*ElouNone))
    
    arrows(leftpteui,ElolMed,leftpteui,ElouMed,length=0.05,angle=90,code=3,lwd=2,lty=3,col="blue")
    #arrows(leftpteui,ElolHigh,leftpteui,ElouHigh,length=0.05,angle=90,code=3,lty=3,lwd=2,col="blue")
    text(leftpteui+0.3,ElolMed-0.02,round(100*ElolMed))
    text(leftpteui+0.3,ElouMed+0.02,round(100*ElouMed))
    #text(leftpteui+0.3,ElolHigh-0.02,round(100*ElolHigh))
    #text(leftpteui+0.3,ElouHigh+0.02,round(100*ElouHigh))
    
    points(rightptii,IhilNone)
    arrows(rightptii,IhilMed,rightptii,IhiuMed,length=0.05,angle=90,code=3,lwd=2)
    #arrows(rightptii,IhilHigh,rightptii,IhiuHigh,length=0.05,angle=90,code=3,lty=2,lwd=2)
    text(rightptii+0.3,IhilNone,round(100*IhilNone))
    text(rightptii+0.3,IhilMed-0.02,round(100*IhilMed))
    text(rightptii+0.3,IhiuMed+0.02,round(100*IhiuMed))
    #text(rightptii+0.3,IhilHigh-0.02,round(100*IhilHigh))
    #text(rightptii+0.3,IhiuHigh+0.02,round(100*IhiuHigh))
    
    arrows(rightptci,EhilNone,rightptci,EhiuNone,length=0.05,angle=90,code=3,lwd=2,lty=2,col="brown")
    text(rightptci+0.3,EhilNone-0.02,round(100*EhilNone))
    text(rightptci+0.3,EhiuNone+0.02,round(100*EhiuNone))
    
    arrows(rightpteui,EhilMed,rightpteui,EhiuMed,length=0.05,angle=90,code=3,lwd=2,lty=3,col="blue")
    #arrows(rightpteui,EhilHigh,rightpteui,EhiuHigh,length=0.05,angle=90,code=3,lty=3,lwd=2,col="blue")
    text(rightpteui+0.3,EhilMed-0.02,round(100*EhilMed))
    text(rightpteui+0.3,EhiuMed+0.02,round(100*EhiuMed))
    #text(rightpteui+0.3,EhilHigh-0.02,round(100*EhilHigh))
    #text(rightpteui+0.3,EhiuHigh+0.02,round(100*EhiuHigh))
    
    legend("bottomright",legend=c("Ign. Int.","95% CI","95% EUI"),
    lty=c(1,2,3),col=c("black","brown","blue"),lwd=3,cex=1.3)
    
    mtext(paste("Point estimate of RRratio = [1 - VE(Low)]/[1 - VE(High)] = ",round(RatioIlNone,2)),side=1,line=5)
    mtext(paste("Ignorance interval for RRratio = ",round(RatioIlMed,2),"--",round(RatioIuMed,2)),side=1,line=6.5)
    mtext(paste("95% CI for RRratio = ",round(RatioElNone,2),"--",round(RatioEuNone,2)),side=1,line=8)
    mtext(paste("95% EUI for RRratio= ",round(RatioElMed,2),"--",round(RatioEuMed,2)),side=1,line=9.5)
    dev.off()
    
}

# UPDATE
# Note: These values are transformed to IU
lowerlimits <- rep(lowerlimitsassays,length(analysispercs)) 
#lowerlimits <- c(0.3076,1.593648,0.242,1.502,
#                 0.3076,1.593648,0.242,1.502,
#                 0.3076,1.593648,0.242,1.502)
markers <- rep(markersforanalysis, length(analysispercs))
# END UPDATE

assay_labels_short
figuretitlesmarkers <- paste0("VE by ",axislabelsmarkers, " > vs. <= ",unlist(cutpoints))
#figuretitlesmarkers      <- c(paste("VE by Day 57 IgG Spike > vs. <= ",cutpoints[1,1]," IU/ml"),
#                              paste("VE by Day 57 IgG RBD   > vs. <= ",cutpoints[2,1],  " IU/ml"),
#                              paste("VE by Day 57 nAb ID50  > vs. <= ",cutpoints[3,1], " IU50/ml"),
#                              paste("VE by Day 57 nAb ID80  > vs. <= ",cutpoints[4,1], " IU80/ml"),
#                              paste("VE by Day 57 IgG Spike > vs. <= ",cutpoints[1,2]," IU/ml"),
#                              paste("VE by Day 57 IgG RBD > vs. <= ",cutpoints[2,2], " IU/ml"),
#                              paste("VE by Day 57 nAb ID50 > vs. <= ",cutpoints[3,2]," IU50/ml"),
#                              paste("VE by Day 57 nAb ID80 > vs. <= ",cutpoints[4,2]," IU80/ml"),
#                              paste("VE by Day 57 IgG Spike > vs. <= ",cutpoints[1,3]," IU/ml"),
#                              paste("VE by Day 57 IgG RBD > vs. <= ",cutpoints[2,3]," IU/ml"),
#                              paste("VE by Day 57 nAb ID50 > vs. <= ",cutpoints[3,3]," IU50/ml"),
#                              paste("VE by Day 57 nAb ID80 > vs. <= ",cutpoints[4,3]," IU80/ml"))

# UPDATE
axislabelsmarkers <- rep(axislabelsmarkers, length(analysispercs))

# This change takes away the colnames of the vector cutpointsmarkers
cutpointsmarkers <- as.vector(cutpoints)
#cutpointsmarkers <- c(cutpoints[1,1],cutpoints[2,1],cutpoints[3,1],cutpoints[4,1],
#                      cutpoints[1,2],cutpoints[2,2],cutpoints[3,2],cutpoints[4,2],
#                      cutpoints[1,3],cutpoints[2,3],cutpoints[3,3],cutpoints[4,3])
# END UPDATE

xlimsset <- list(NA,length(axislabelsmarkers))
ataxisset <- list(NA,length(axislabelsmarkers))
labelsxaxisset <- list(NA,length(axislabelsmarkers))

for (i in 1:length(axislabelsmarkers)) {
    xlimsset[[i]] <- c(-1,6)
    ataxisset[[i]] <- c(log10(lowerlimits[i]),1,2,3,4,5)
    labelsxaxisset[[i]] <- c("LOD",expression(10),expression(10^2),
    expression(10^3),expression(10^4),expression(10^5))
}


# UPDATE
filenamesplots <- paste0(save.results.to, "plotPSbin", study_name, "Day", tpeak, assays, "cut", rep(1:length(analysispercs), each=length(assays)),".pdf")
# END UPDATE

#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57Spikecut1.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57RBDcut1.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57ID50cut1.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57ID80cut1.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57Spikecut2.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57RBDcut2.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57ID50cut2.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57ID80cut2.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57Spikecut3.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57RBDcut3.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57ID50cut3.pdf",
#"T:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/reports/plots/plotPSbinCOVEDay57ID80cut3.pdf")

#filenamesplots <- c("H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57Spikecut1.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57RBDcut1.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57ID50cut1.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57ID80cut1.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57Spikecut2.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57RBDcut2.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57ID50cut2.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57ID80cut2.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57Spikecut3.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57RBDcut3.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57ID50cut3.pdf",
#"H:/Coronavirus/HVTNstudies/ModernaProgram/ImmuneCorrelatesManuscript/PrincipalStratificationBinary/plotPSbinCOVEDay57ID80cut3.pdf")

#filenamesplots <- c("C:/temp/plotPSbinCOVEDay57Spikecut1.pdf",
#"C:/temp/plotPSbinCOVEDay57RBDcut1.pdf",
#"C:/temp/plotPSbinCOVEDay57ID50cut1.pdf",
#"C:/temp/plotPSbinCOVEDay57ID80cut1.pdf",
#"C:/temp/plotPSbinCOVEDay57Spikecut2.pdf",
#"C:/temp/plotPSbinCOVEDay57RBDcut2.pdf",
#"C:/temp/plotPSbinCOVEDay57ID50cut2.pdf",
#"C:/temp/plotPSbinCOVEDay57ID80cut2.pdf",
#"C:/temp/plotPSbinCOVEDay57Spikecut3.pdf",
#"C:/temp/plotPSbinCOVEDay57RBDcut3.pdf",
#"C:/temp/plotPSbinCOVEDay57ID50cut3.pdf",
#"C:/temp/plotPSbinCOVEDay57ID80cut3.pdf")


# Make the plots for the set of markers:
for (i in 1:length(labelsmarkersforanalysis)) {
    
    markernoncases <- data[Delta.Day1==0,markers[i]]
    markername <- labelsmarkersforanalysis[i]
    xlab1 <- axislabelsmarkers[i]
    filenameplot <- filenamesplots[i] 
    main1 <- figuretitlesmarkers[i]
    lowerlimit <- lowerlimits[i] 
    cutpt <- log10(unlist(cutpointsmarkers)[i])
    xlims <- xlimsset[[i]]
    atxaxis <- ataxisset[[i]]
    labelsxaxis <- labelsxaxisset[[i]]
    
    leftptii <- -1.0
    leftptci <- -0.3
    leftpteui <- 0.4
    rightptii <- 3.6
    rightptci <- 4.3
    rightpteui <- 5.0
    plotprincstratbinary(datareport,markernoncases,markername,xlims,xlab1,filenameplot,main1,atxaxis,labelsxaxis,cutpt,leftptii,leftptci,leftpteui,rightptii,rightptci,rightpteui) 
    
}
