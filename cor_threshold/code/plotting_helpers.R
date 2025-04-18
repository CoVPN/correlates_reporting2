
#' Generates plot of threshold-response function.
#' @param marker The marker variable to generate plots for.
#' @param simultaneous_CI True if simultaneous CI should be plotted. Otherwise if False pointwise CI are plotted.
#' @param monotone True if monotone correction should be done on estimates. Otherwise no monotone correction. This should be done if one expects monotonicity.
#' Assume monotone nonincreasing
get_plot <- function(marker, simultaneous_CI = F, monotone = T, above = TRUE) {
  
  risk_plac = prev.plac # computed in _common.R
  
  is.delta=startsWith(marker,"Delta")
  
  key <- marker
  if(above){
    append <- ""
  } else {
    append <- "_below"
  }
  if(monotone) {
    load(file = here::here("output", TRIAL, COR, paste0("tmleThresh_monotone_", marker,append, ".RData")))
  } else {
    load(file = here::here("output", TRIAL, COR, paste0("tmleThresh_",  marker,append, ".RData")))
  }
  str(esttmle)
  
  
  time <- tpeak
  day <- ""
  if(TRIAL == "hvtn705second"){
    laby <- paste0("Probability of HIV by Day ",max_t)
  } else {
    laby <- paste0("Probability of COVID by Day ",max_t)
  }
  
  labx <- plotting_assay_label_generator(marker, above)
  # subtitle_main <- "Nonparametric estimate of Threshold-response function"
  data <- read.csv(here::here("output", TRIAL, COR, "data_clean", paste0("data_secondstage.csv")))
  main <- plotting_assay_title_generator(marker)
  if(length(grep("start", key)) > 0) {
    main <- paste0(main , " (1-day-post)")
  }
  print(main)
  #paste0("Cumulative Risk of COVID by Day ", tf[time])
  
  ident <- function(x) x
  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)
  # hist(na.omit(data_biased[[marker]]),col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  # Get initial threshold-response plot with simultaneous CI
  v <- plot_threshold_response(esttmle, simultaneous_CI = simultaneous_CI, monotone = F)
  data_tmp <- na.omit(data[, c(marker, "Ttilde", "Delta", "wt"), drop = F])
  max_thresh <- max(data_tmp[data_tmp[["Delta"]] == 1, marker])
  # out <- hist(na.omit(data_tmp[[marker]]), col=col,axes=F,labels=F,main="",xlab="",ylab="",breaks=10,border=0,freq=F)
  scale_coef <- max(v$data$upper, na.rm = T) * 1
  
  # coef <- max(out$density)/scale_coef
  RCDF <- function(a) {
    sum(data_tmp$wt * (data_tmp[[marker]] >= a)) / sum(data_tmp$wt) * scale_coef
  }
  RCDF <- Vectorize(RCDF)
  if (!simultaneous_CI) {
    #main <- paste0(main, " with point-wise confidence intervals")
  } else {
    #main <- paste0(main, " with simultaneous confidence bands")
  }
  a <- marker_to_assay[[marker]]
  print(marker)
  print(a)
  if (!is.delta) {
    xlim=get.range.cor(data, a, tpeak)
    
    # hack
    if (TRIAL=="moderna_boost" & COR=="BD29nnaive") {
      # use naive data to set xlim
      data1 <- read.csv(here::here("output", TRIAL, "BD29naive", "data_clean", paste0("data_secondstage.csv")))
      xlim=get.range.cor(data1, a, tpeak)
    }
    
  } else {
    xlim=range(data[[marker]], na.rm=T)
    
    # hack
    if (TRIAL=="moderna_boost" & COR=="BD29nnaive") {
      # use naive data to set xlim
      data1 <- read.csv(here::here("output", TRIAL, "BD29naive", "data_clean", paste0("data_secondstage.csv")))
      xlim=range(data1[[marker]], na.rm=T)
    }
  }
  
  
  print(quantile(data[[marker]]))
  print(xlim)
  llod <- lloxs[a]
  labels_info <- draw.x.axis.cor(xlim, lloxs[a], if(is.delta) "delta" else llox_labels[a], for.ggplot=T)
  print(labels_info)
  xx <- labels_info$ticks
  labels <- as.list(labels_info$labels)
  
  xlimits <- xlim
  print(xlimits)
  
  print(main)
  print(v) 
  
  plot <- v + ggtitle(main) +
    stat_function(fun = RCDF, color = col, geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      name = laby,
      sec.axis = sec_axis(~ . / scale_coef, name = "Reverse CDF"), n.breaks = 10
    )  +
    theme(plot.title = element_text(size = 25), axis.text.x = element_text(angle = 0, hjust = 1, size = 18), axis.text.y = element_text(angle = 0, hjust = 1, size = 18)) +
    # geom_hline(aes(yintercept=risk_vac), alpha = 0.4) + geom_text(alpha = 0.75,aes(median(v$data$cutoffs),risk_vac,label = "vaccine overall risk"), vjust = -0.5, size = 5) +
    # comment out till we get risk_plac
    # geom_text(alpha = 0.75, aes(quantile(v$data$cutoffs, 0.1),min(max(v$data$upper),risk_plac),label = paste0("placebo overall risk: ", risk_plac)), vjust = 0, size = 5) +
    scale_x_continuous(
      breaks = xx,#union(floor(esttmle[, 1]), ceiling(esttmle[, 1])),
      labels = do.call(expression,labels),
      name =labx,
      limits = xlimits
      #trans_format("ident", math_format(10^.x)),
      #limits = c(min(esttmle[, 1]) - 0.1, max(esttmle[, 1]) + 0.1)
    )
  
  if(above  && max_thresh < log10(uloqs[a]) - 0.05) {
		  # comment out the following line to skip red lines since it shows up in some plots and not others and it is not clear what the x values are for the red lines
#    plot <- plot + geom_vline(xintercept = max_thresh, colour = "red", linetype = "longdash")
  } 
  # comment out till we get risk_plac
  # else if(!above && risk_plac<= max(v$data$upper, na.rm = T)) {
  #   plot <- plot + geom_hline(aes(yintercept=risk_plac), alpha = 0.4, linetype = "longdash", colour = "red")
  # }
  plot
  
  data_interp <- as.data.frame(approx(plot$data$cutoffs, plot$data$est, xout = data[data$Delta==1,marker], rule = 2  ))
  
  
  plot <- plot  +geom_point(data=data_interp,aes(x=x, y=y), colour = "blue")
  
  #+  geom_text(aes(x=max_thresh *(1.01), label="No observed events", y=0.002), colour="black", angle=90, text=element_text(size=11))
  append_end <- ""
  append_start <- "PLOT"
  folder <- ""
  if (monotone) {
    append_start <- paste0(append_start, "_monotone_")
  } else {
    append_start <- paste0(append_start, "_")
  }
  if (simultaneous_CI) {
    append_end <- paste0(append_end, "_", "simultCI",append)
    folder <- "simultaneous_CI"
  } else {
    append_end <- paste0(append_end, "_", "pointwiseCI",append)
    folder <- "pointwise_CI"
  }
  
  print(folder)
  
  ggsave(
    filename = here::here("output", TRIAL, COR,
      "figs", folder,
      paste0(append_start, marker, append_end, ".pdf")
    ),
    plot = plot, height = 7, width = 9
  )
  
  return(plot)
}

# Generates tables (both pointwise CI and simultaneous CI) for Threshold-response function estimates
#' @param marker The marker variable to generate plots for.
#' @param num_show The number of thresholds to include in table.
generate_tables <- function(marker, num_show = 10, monotone = F, above = T) {
  print("TABLES")
  
  data_secondstage <- read.csv(here::here("output", TRIAL, COR, "data_clean", paste0("data_secondstage.csv")))
  
  if(above){
    append <- ""
  } else {
    append <- "_below"
  }
  if(monotone) {
    load(file = here::here("output", TRIAL, COR, paste0("tmleThresh_monotone_",  marker,append, ".RData")))
  } else {
    load(file = here::here("output", TRIAL, COR, paste0("tmleThresh_",  marker,append, ".RData")))
  }
  esttmle_table <- esttmle
  esttmle_table[, 1] <- round(esttmle_table[, 1], 3)
  esttmle_table[, 2] <- round(esttmle_table[, 2], 5)
  esttmle_table[, 3] <- round(esttmle_table[, 3], 5)
  esttmle_table[, 4] <- round(esttmle_table[, 4], 5)
  esttmle_table[, 5] <- round(esttmle_table[, 5], 5)
  esttmle_table[, 6] <- round(esttmle_table[, 6], 5)
  esttmle_table[, 7] <- round(esttmle_table[, 7], 5)
  esttmle_table <- esttmle_table[,c(1:7)]
  esttmle_table <- data.frame(esttmle_table, paste0(gsub("e\\+0|e\\-0", " * 10$^{", format(10^esttmle_table[, 1], scientific = T, digits = 3)), "}$"))
  esttmle_table <- esttmle_table[, c(1, 8, 2, 4, 5, 6, 7)]
  esttmle_table[esttmle_table[, 3] < 0, 3] <- 0
  colnames(esttmle_table) <- c("log$_{10}$-Threshold", "Threshold", "Risk estimate", "CI left", "CI right", "CI left", "CI right")
  # Save nice latex table
  thresh_mand <- report.assay.values(data_secondstage[[marker]], marker)
  index_to_show <- sapply(thresh_mand, function(thresh) {
    which.min(abs(thresh-esttmle_table[,1]))
  })
  #index_to_show <- unique(round(seq.int(1, nrow(esttmle_table), length.out = num_show)))
  
  ptwise_tab_guts <- esttmle_table[index_to_show, c(1, 2, 3, 4, 5)]
  
  if(monotone) {
    aadd1 <- paste0(append,"monotone_")
  } else {
    aadd1 <- append
  }
  saveRDS(ptwise_tab_guts, file = here::here("output", TRIAL, COR,
    "figs", "pointwise_CI",
    paste0("TABLE_",aadd1,  marker, "_pointwiseCI.rds")
  ))
  simul_tab_guts <- esttmle_table[index_to_show, c(1, 2, 3, 6, 7)]
  
  saveRDS(simul_tab_guts, file = here::here("output", TRIAL, COR,
    "figs", "simultaneous_CI",
    paste0("TABLE_", aadd1, marker, "_simultCI.rds")
  ))
  return(list(pointwise = esttmle_table[index_to_show, c(1, 2, 3, 4, 5)], simult = esttmle_table[index_to_show, c(1, 2, 3, 6, 7)]))
} 
