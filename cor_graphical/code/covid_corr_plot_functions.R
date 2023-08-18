#' a sub-function called by function: ggally_statistic_resample
#' 
#' when B > 1, resamping-based correlation, used in ggplots, allowing for strata (if no strata, need to input all 1's as strata)
#' when B = 0, weighted correlation, used in ggplots, allowing for strata (if no strata, need to input all 1's as strata)
ggally_statistic_resample <- function(
  data,
  mapping,
  B = 0, # when B = 0, no resampling will be done, and weighted spearman correlation will be done
  strata = NULL,
  weight = NULL,
  text_fn,
  title,
  na.rm = NA,
  display_grid = FALSE,
  justify_labels = "right",
  justify_text = "left",
  sep = ": ",
  family = "mono",
  title_args = list(),
  group_args = list(),
  align_percent = 0.5,
  title_hjust = 0.5,
  group_hjust = 0.5) {
  set_if_not_there <- function(obj, key, value) {
    obj <- as.list(obj)
    # if (! "family" %in% rlang::names2(obj)) {
    #  obj$family <- family
    # }
    obj
  }
  
  # title_args <- set_if_not_there(title_args, "family", family)
  # group_args <- set_if_not_there(group_args, "family", family)
  
  title_args <- set_if_not_there(title_args, "hjust", title_hjust)
  group_args <- set_if_not_there(group_args, "hjust", group_hjust)
  
  xData <- eval_data_col(data, mapping$x)
  yData <- eval_data_col(data, mapping$y)
  strataData <- strata
  weightData <- weight
  colorData <- eval_data_col(data, mapping$colour)
  
  if (is.numeric(colorData)) {
    stop("`mapping` color column must be categorical, not numeric")
  }
  
  display_na_rm <- is.na(na.rm)
  if (display_na_rm) {
    na.rm <- TRUE
  }
  if (isTRUE(na.rm)) {
    if (!is.null(colorData) && (length(colorData) == length(xData))) {
      rows <- complete.cases(xData, yData, colorData, strataData, weightData)
    } else {
      rows <- complete.cases(xData, yData, strataData, weightData)
    }
    
    if (any(!rows)) {
      if (!is.null(colorData) && (length(colorData) == length(xData))) {
        colorData <- colorData[rows]
      }
      xData <- xData[rows]
      yData <- yData[rows]
      strataData <- strataData[rows]
      weightData <- weightData[rows]
      
      if (isTRUE(display_na_rm)) {
        total <- sum(!rows)
        if (total > 1) {
          warning("Removed ", total, " rows containing missing values")
        } else if (total == 1) {
          warning("Removing 1 row that contained a missing value")
        }
      }
    }
  }
  xVal <- xData
  yVal <- yData
  strataVal <- strataData
  weightVal <- weightData
  # if the mapping has to deal with the data, remove it
  ### NOTE: IDK what this does. inherited from old code.
  for (mappingName in names(mapping)) {
    itemData <- eval_data_col(data, mapping[[mappingName]])
    if (!inherits(itemData, "AsIs")) {
      mapping[[mappingName]] <- NULL
    }
  }
  ### END IDK
  
  # calculate variable ranges so the gridlines line up
  xValNum <- as.numeric(xVal)
  yValNum <- as.numeric(yVal)
  
  xmin <- min(xValNum, na.rm = TRUE)
  xmax <- max(xValNum, na.rm = TRUE)
  xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 * (xmax - xmin))
  ymin <- min(yValNum, na.rm = TRUE)
  ymax <- max(yValNum, na.rm = TRUE)
  yrange <- c(ymin - 0.01 * (ymax - ymin), ymax + 0.01 * (ymax - ymin))
  # if there is a color grouping...
  if (!is.null(colorData) && !inherits(colorData, "AsIs")) {
    cord <- ddply(
      data.frame(
        x = xData, y = yData, st = strataData, wt = weightData,
        color = colorData
      ),
      "color",
      function(dt) {
        text_fn(dt$x, dt$y, dt$st, dt$wt, B)
      }
    )
    colnames(cord)[2] <- "text"
    
    # put in correct order
    lev <- levels(as.factor(colorData))
    ord <- rep(-1, nrow(cord))
    for (i in 1:nrow(cord)) {
      for (j in seq_along(lev)) {
        if (identical(as.character(cord$color[i]), as.character(lev[j]))) {
          ord[i] <- j
        }
      }
    }
    cord <- cord[order(ord[ord >= 0]), ]
    
    # make labels align together
    cord$label <- str_c(
      format(cord$color, justify = justify_labels),
      sep,
      format(cord$text, justify = justify_text)
    )
    
    # title
    ggally_text_args <- append(
      list(
        label = str_c(title, sep, text_fn(xVal, yVal, strataVal, wVal, B)),
        mapping = mapping,
        xP = 0.5,
        yP = 0.9,
        xrange = xrange,
        yrange = yrange
      ),
      title_args
    )
    p <- do.call(ggally_text, ggally_text_args)
    
    xPos <- rep(align_percent, nrow(cord)) * diff(xrange) +
      min(xrange, na.rm = TRUE)
    yPos <- seq(
      from = 0.9,
      to = 0.2,
      length.out = nrow(cord) + 1
    )
    yPos <- yPos * diff(yrange) + min(yrange, na.rm = TRUE)
    yPos <- yPos[-1]
    
    cordf <- data.frame(xPos = xPos, yPos = yPos, labelp = cord$label)
    cordf$labelp <- factor(cordf$labelp, levels = cordf$labelp)
    
    # group text values
    geom_text_args <- append(
      list(
        data = cordf,
        aes(
          x = xPos,
          y = yPos,
          label = labelp,
          color = labelp
        )
      ),
      group_args
    )
    p <- p + do.call(geom_text, geom_text_args)
  } else {
    ggally_text_args <- append(
      list(
        label = paste0(title, sep, text_fn(
          xVal, yVal, strataVal,
          weightVal, B
        ), collapse = ""),
        mapping,
        xP = 0.5,
        yP = 0.5,
        xrange = xrange,
        yrange = yrange
      ),
      title_args
    )
    p <- do.call(ggally_text, ggally_text_args)
  }
  
  if (!isTRUE(display_grid)) {
    p <- p +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
          linetype = "solid",
          color = theme_get()$panel.background$fill,
          fill = "transparent"
        )
      )
  }
  p + theme(legend.position = "none")
}


#' a sub-function called by function: covid_corr_pairplots
#' 
#' when B > 1, resamping-based correlation, used in ggplots, allowing for strata (if no strata, need to input all 1's as strata)
#' when B = 0, weighted correlation, used in ggplots, allowing for strata (if no strata, need to input all 1's as strata)
ggally_cor_resample <- function(
  data,
  mapping,
  strata,
  weight,
  B = 0, # if B == 0, then no resampling will be done and weighted spearman correlation coef will be computed
  # if B > 1, then resamping will be done B times and simple spearman correlated coef will be computed and averaged
  seed = 12345,
  ...,
  stars = TRUE,
  method = "spearman",
  use = "complete.obs",
  display_grid = FALSE,
  digits = 3,
  title_args = list(...),
  group_args = list(...),
  justify_labels = "right",
  align_percent = 0.5,
  title = "Corr",
  alignPercent = warning("deprecated. Use `align_percent`"),
  displayGrid = warning("deprecated. Use `display_grid`")) {
  if (!missing(alignPercent)) {
    warning("`alignPercent` is deprecated. Please use `align_percent` if alignment still needs to be adjusted")
    align_percent <- alignPercent
  }
  if (!missing(displayGrid)) {
    warning("`displayGrid` is deprecated. Please use `display_grid`")
    display_grid <- displayGrid
  }
  
  na.rm <-
    if (missing(use)) {
      # display warnings
      NA
    } else {
      (use %in% c("complete.obs", "pairwise.complete.obs", "na.or.complete"))
    }
  
  ggally_statistic_resample(
    data = data,
    mapping = mapping,
    strata = strata,
    weight = weight,
    na.rm = na.rm,
    align_percent = align_percent,
    display_grid = display_grid,
    title_args = title_args,
    group_args = group_args,
    justify_labels = justify_labels,
    justify_text = "left",
    sep = if ("colour" %in% names(mapping)) ": " else ":\n",
    title = title,
    text_fn = function(x, y, st, wt, B) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      nn <- length(x)
      corvec <- rep(NA, B)
      set.seed(seed)
      resamp_mat <- sapply(1:B, function(ii) sample.int(n = nn, replace = TRUE, prob = wt))
      # write.csv(data.frame(x = x, y = y, strata = st), "input_columns.csv", row.names = FALSE)
      
      # write.csv(resamp_mat, "output_row_number.csv", row.names = FALSE)
      if (B > 1) { # if B > 1, resampling will be done
        for (bb in seq_len(B)) {
          resamp_vec <- resamp_mat[, bb]
          x_resamp <- x[resamp_vec]
          y_resamp <- y[resamp_vec]
          st_resamp <- st[resamp_vec]
          
          # write.csv(data.frame(x = x_resamp, y = y_resamp, strata = st_resamp), "input_columns_last_resamp.csv", row.names = FALSE)
          
          suppressWarnings(st_resamp_dummy <-
                             dummies::dummy(st_resamp, sep = "_"))
          
          resamp_data <- cbind(
            data.frame(x = x_resamp, y = y_resamp),
            st_resamp_dummy[, 1:(ncol(st_resamp_dummy) - 1)]
          )
          names(resamp_data)[3:ncol(resamp_data)] <-
            paste0("strata", 1:(ncol(st_resamp_dummy) - 1))
          #write.csv(resamp_data, "resamp_data.csv", row.names = FALSE)
          
          #fml <- formula(paste0(
          #    "y | x ~ ",
          #    paste0("strata", 1:(ncol(st_resamp_dummy) - 1),
          #           collapse = "+"
          #    )
          #))
          
          #if (study_name=="COVEBoost"){ 
          # e.g. PROFISCOV's Bstratum variable only includes one value, so use simple spearman for this study, need to specify "exact = FALSE" because of ties
          corObj <- try(cor.test(resamp_data$x, resamp_data$y, method = "spearman", exact = FALSE)$estimate,
                        silent = TRUE
          )
          
        }
        #} #else { # Partial Spearman's Rank Correlation
        #  corObj <- try(partial_Spearman(
        #      formula = fml, data = resamp_data,
        #      fit.x = "lm", fit.y = "lm"
        #  )$TS$TB$ts,
        #  silent = TRUE
        #  )
        #}
        
        # make sure all values have X-many decimal places
        corvec[bb] <- ifelse(class(corObj) == "try-error",
                             NA,
                             corObj
        )
      } else if (B == 0){ # resampling will not be done, weighted spearman correlation will be done
        
        corObj <- try(weightedCorr(x, y, method = "Spearman", weights = wt),
                      silent = TRUE
        )
        
        corvec <- ifelse(class(corObj) == "try-error",
                         NA,
                         corObj
        )
        
      }
      #saveRDS(corvec, file = "corvec.RDS")
      cor_est <- mean(corvec, na.rm = TRUE)
      cor_txt <- formatC(cor_est, digits = digits, format = "f")
      cor_txt
    }
  )
}

#' Pairplots of assay readouts
#'
#' Produce the pairplots of assay readouts. The correlation is calculated by
#' the resampling-based strata adjusted Spearman rank correlation
#' @param plot_dat: data frame: data for plotting.
#' @param time: string: one of "D1", "D29", "D57", "Delta29overB" or
#'  "Delta57overB".
#' @param assays: vector of strings: the assay names for plotting.
#' @param strata: string: the column name in plot_dat that indicates the
#'  strata. when no strata is needed, a column of 1's needs to be assigned
#' @param weight: string: the column name in plot_dat that indicates the
#'  individual sampling weights.
#' @param plot_title: string: title of the plot.
#' @param column_labels: vector of strings: titles of each column.
#' @param height: scalar: plot height.
#' @param width: scalar: plot width.
#' @param units: string: the unit of plot height and width.
#' @param corr_size: scalar: font size of the correlation labels.
#' @param point_size: scalar: point size in the scatter plots.
#' @param loess_lwd: scalar: loess line width in the scatter plots.
#' @param plot_title_size: scalar: font size of the plot title.
#' @param column_label_size: scalar: font size of the column labels.
#' @param axis_label_size: scalar: font size of the axis labels.
#' @param filename: string: output file name.
#' @param write_to_file: logical: whether to output file or just output an object
#'
#' @return pairplots: a ggplot object of the pairplot
covid_corr_pairplots <- function(plot_dat, ## data for plotting
                                 time,
                                 assays,
                                 strata, 
                                 weight,
                                 plot_title,
                                 column_labels,
                                 seed = 12345,
                                 height = 6.1,
                                 width = 6.05,
                                 units = "in",
                                 corr_size = 5,
                                 point_size = 0.5,
                                 loess_lwd = 1,
                                 plot_title_size = 10,
                                 column_label_size = 6.5,
                                 axis_label_size = 9,
                                 filename,
                                 write_to_file = T) {
  dat.tmp <- plot_dat[, paste0(time, assays)]
  rr <- range(dat.tmp, na.rm = TRUE)
  
  if (rr[1] == rr[2]) {
    rr <- c(rr[1] - 1, rr[2] + 1)
  }
  
  if (rr[2] - rr[1] < 2) {
    rr <- floor(rr[1]):ceiling(rr[2])
  }
  
  breaks <- floor(rr[1]):ceiling(rr[2])
  
  # readjust limits to include all of the breaks
  rr <- c(floor(rr[1]), ceiling(rr[2]))
  
  if (max(breaks) - min(breaks) >= 6) {
    breaks <- breaks[breaks %% 2 == 0]
  }
  
  pairplots <- ggpairs(
    data = dat.tmp, title = plot_title,
    columnLabels = column_labels,
    upper = list(
      continuous =
        wrap(ggally_cor_resample,
             stars = FALSE,
             size = corr_size,
             seed = seed,
             strata = plot_dat[, strata],
             weight = plot_dat[, weight]
        )
    ),
    lower = list(
      continuous =
        wrap("points", size = point_size)
    )
  ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = plot_title_size),
      strip.text = element_text(size = column_label_size, face = "bold"),
      strip.background = element_rect(fill=NA,colour=NA),
      axis.text = element_text(size = axis_label_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
  pairplots[1, 1] <- pairplots[1, 1] +
    scale_x_continuous(limits = rr, breaks = breaks) + ylim(0, 1.2)
  for (j in 2:pairplots$nrow) {
    for (k in 1:(j - 1)) {
      pairplots[j, k] <- pairplots[j, k] +
        stat_smooth(
          method = "loess", color = "red", se = FALSE,
          lwd = loess_lwd
        ) +
        scale_x_continuous(
          limits = rr, breaks = breaks,
          labels = scales::math_format(10^.x)
        ) +
        scale_y_continuous(
          limits = rr, breaks = breaks,
          labels = scales::math_format(10^.x)
        )
    }
    pairplots[j, j] <- pairplots[j, j] +
      scale_x_continuous(
        limits = rr, breaks = breaks,
        labels = scales::math_format(10^.x)
      ) + ylim(0, 1.2)
  }
  
  if (write_to_file == T){
    ggsave(
      filename = filename, plot = pairplots, width = width, height = height,
      units = units
    )
  } else{
    pairplots
  }
}
