#-----------------------------------------------
# obligatory to append to the top of each script
#-----------------------------------------------

#' Generate response calls for endpoints:
#' Responders for ID50/80 neutralization titer at each pre-defined timepoint are defined as participants who
#' had baseline values below the LLOD with detectable ID50/80 neutralization titer
#' above the assay LLOD, or as participants with baseline values above
#' the LLOD with a 4-fold increase in neutralizing antibody titer.
#'
#' Responders for bAb at each pre-defined timepoint are defined as participants who had values above a certain value, i.e., pos.cutoff.
#' 
#' @param data The dataframe with the endpoint of interest at baseline and
#'  post-baseline.
#' @param post_times The variable of endpoint value at post-baseline
#' @param assays Assay list that needs a response call
#' @param respoderFR Fold-rise used for response call positivity criteria
#' @param pos.cutoffs List of positivity call thresholds for all assays in the assay list
#'
#' @return Response calls for \code{bl} and \code{post}
getResponder <- function(data,
                         post_times, 
                         assays,
                         responderFR = 4,
                         pos.cutoffs = pos.cutoffs) {
    
    for (i in post_times){
        for (j in assays){
            post <- paste0("Day", i, j)
            bl <- paste0("B", j)
            delta <- paste0("Delta", i, "overB", j)
            
            if (grepl("bind", j)) {
                
            data[, paste0(post, "Resp")] <- as.numeric(data[, post] > log10(pos.cutoffs[j]))
            data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(pos.cutoffs[j]))
            
            } else { 
                data[, paste0(post, "Resp")] <- as.numeric(
                    (data[, bl] < log10(pos.cutoffs[j]) & data[, post] > log10(pos.cutoffs[j])) |
                        (data[, bl] >= log10(pos.cutoffs[j]) & as.numeric(10^data[, delta] >= responderFR)))
                data[, paste0(bl, "Resp")] <- as.numeric(data[, bl] > log10(pos.cutoffs[j]))
            }
            
        }
    }
    return(data)
}


#' response rate is calculated using weight
#'
#'
#' @param data The dataframe in the long format of assays and times, with response variable
#' @param group The group level at which response rate should be calculated
#'
#' @return summary statistics and response rate by group
# 
get_desc_by_group <- function(data, 
                              group=group){
    
    data$wt = data$wt.subcohort
    
    complete <- complete.cases(data[, group])
    
    dat_stats <-
        data %>% filter(complete==1) %>%
        group_by_at(group) %>%
        mutate(counts = n(),
               num = sum(response * wt, na.rm=T),
               denom = sum(wt, na.rm=T),
               #N_RespRate = paste0(counts, "\n",round(num/denom*100, 1),"%"),
               RespRate = ifelse(denom!=0 && !is.na(pos.cutoffs), paste0(round(num/denom*100, 1),"%"), ""), # RespRate at Delta timepoints will be ""
               min = min(value),
               q1 = quantile(value, 0.25, na.rm=T),
               median = median(value, na.rm=T),
               q3 = quantile(value, 0.75, na.rm=T),
               max= max(value))
    
    return(dat_stats)
}

#' A ggplot object for violin box plot without lines, loop by timepoints
#' 
#' @param dat Dataframe in long format of assays and timepoints
#' @param assays List of assays for plots
#' @param times List of times for plots
#' @param ylim y-axis limit
#' @param panel.text.size font size for text within panels
#' @param axis.x.text.size font size for x-axis tick label
#' @param strip.x.text.size font size for x-axis strip label
#' @param facet.x.var horizontal facet variable 
#' @param facet.y.var vertical facet variable 
#' @return A ggplot object list for violin + box plot without lines
f_by_time_assay <- 
    function(dat,
             assays = assays,
             times = times,
             ylim = c(0,7.2),
             panel.text.size = 2.2,
             axis.x.text.size = 18,
             strip.x.text.size = 10,
             facet.y.var,
             facet.x.var
    ) {
        
        plot_theme <- theme_bw(base_size = 25) +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size = axis.x.text.size),
                  axis.text.y = element_text(size = 25),
                  axis.title = element_text(size = 24, face="bold"),
                  strip.text.x = element_text(size = strip.x.text.size), # facet label size
                  strip.text.y = element_text(size = 25),
                  axis.ticks.x=element_blank(),
                  strip.background = element_rect(fill=NA,colour=NA),
                  strip.placement = "outside",
                  legend.position = "bottom", 
                  legend.text = element_text(size = 16, face="plain"),
                  legend.key = element_blank(), # remove square outside legend key
                  plot.caption = element_text(size = 26, hjust=0, face="plain"), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.margin = margin(5.5, 12, 5.5, 5.5, "pt")) 
        
        p1 <- dat %>%
            filter(assay %in% assays & time %in% times) %>%
            left_join(assay_metadata, by="assay") %>%
            mutate(Trt_nnaive = factor(paste(Trt, nnaive), 
                                       levels = c("Vaccine Naive", "Vaccine Non-naive", "Placebo Naive", "Placebo Non-naive"),
                                       labels = c("Vaccine\nnaive", "Vaccine\nnon-naive", "Placebo\nnaive", "Placebo\nnon-naive")),
                   assay_label_short = gsub("PsV Neutralization to |PsV Neutralization |Binding Antibody to Spike ", "", assay_label)
            ) %>%
            ungroup() %>%
            group_split(time) %>%
            purrr::map(function(d){
                ggplot(data = d, aes(x = x, y = value)) +
                    facet_grid(rows = facet.y.var, cols = facet.x.var) +
                    geom_violin(scale = "width", na.rm = TRUE, show.legend = FALSE, color = "#FF6F1B") +
                    geom_boxplot(width = 0.25, lwd = 1.5, alpha = 0.3, stat = "boxplot", outlier.shape = NA, show.legend = FALSE, color = "#FF6F1B") +
                    geom_jitter(width = 0.1, height = 0, size = 2, show.legend = TRUE, color = "#FF6F1B") +
                    #scale_color_manual(name = "", values = "#FF6F1B", guide = "none") + # guide = "none" in scale_..._...() to suppress legend
                    # The lower and upper hinges correspond to the first and third quartiles (the 25th and 75th percentiles)
                    # Whisker: Q3 + 1.5 IQR
                    geom_text(aes(label = ifelse(RespRate!="","Rate",""), x = 0.4, y = 5), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE) +
                    geom_text(aes(x = x, label = RespRate, y = 5), color = "black", size = panel.text.size, check_overlap = TRUE) +
                    
                    geom_hline(aes(yintercept = ifelse(RespRate!="",lbval,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
                    geom_text(aes(label = ifelse(RespRate!="",lb,""), x = 0.4, y = lbval), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE, na.rm = TRUE) + 
                    # only plot uloq for ID50
                    geom_hline(aes(yintercept = ifelse(RespRate!="",lbval2,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
                    geom_text(aes(label = ifelse(RespRate!="",lb2,""), x = 0.4, y = lbval2), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE, na.rm = TRUE) + 
                    scale_x_discrete(labels = "") + 
                    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], 3), labels = scales::math_format(10^.x)) +
                    labs(x = "Assay", y = unique(d$panel), title = paste(unique(d$panel), "distributions at", unique(d$time)), color = "Category", shape = "Category") +
                    plot_theme +
                    guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
            })
        return(p1)
    }


#' A ggplot object for longitudinal violin box plot with lines, loop by assays
#' 
#' @param dat Dataframe in long format of assays and timepoints
#' @param x.var x variable, e.g. time_cohort, time
#' @param x.lb x variable label
#' @param assays List of assays for plots
#' @param times List of times for plots
#' @param ylim y-axis limit
#' @param ybreaks y-axis breaks
#' @param panel.text.size font size for text within panels
#' @param facet.x.var horizontal facet variable 
#' @param facet.y.var vertical facet variable
#' @param strip.text.y.size strip label size for y-axis, default is 25
#' @param strip.text.x.size strip label size for x-axis, default is 25
#' @param axis.text.x.size x-axis label size, default is 9.5
#' @return A ggplot object list for longitudinal violin + box plot with lines
f_longitude_by_assay <- function(
    dat,
    x.var = "time",
    x.lb = c("D1","D22","D43"),
    assays = assays,
    times = times,
    ylim = c(0,7.2),
    panel.text.size = 4,
    facet.x.var,
    facet.y.var,
    strip.text.y.size = 25,
    strip.text.x.size = 25,
    axis.text.x.size = 15
) {
    
    plot_theme <- theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(size = axis.text.x.size),
              axis.text.y = element_text(size = 18),
              axis.title = element_text(size = 24, face="bold"),
              strip.text.x = element_text(size = strip.text.x.size), # facet label size
              strip.text.y = element_text(size = strip.text.y.size),
              strip.background = element_rect(fill=NA,colour=NA),
              strip.placement = "outside",
              legend.position = "bottom", 
              legend.text = element_text(size = 16, face="plain"),
              legend.key = element_blank(), # remove square outside legend key
              plot.caption = element_text(size = 26, hjust=0, face="plain"), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.margin = margin(5.5, 12, 5.5, 5.5, "pt")) 
    
    p2 <- dat %>%
        filter(assay %in% assays & time %in% times) %>%
        left_join(assay_metadata, by="assay") %>%
        mutate(Trt_nnaive = factor(paste(Trt, nnaive), 
                                   levels = c("Vaccine Naive", "Vaccine Non-naive", "Placebo Naive", "Placebo Non-naive"),
                                   labels = c("Vaccine\nnaive", "Vaccine\nnon-naive", "Placebo\nnaive", "Placebo\nnon-naive")),
               assay_label_short = gsub("PsV Neutralization to |PsV Neutralization |Binding Antibody to Spike ", "", assay_label)
        ) %>%
        ungroup() %>%
            ggplot(aes_string(x = x.var, y = "value")) +
                facet_grid(rows = facet.y.var, col = facet.x.var) +
                
                geom_violin(scale = "width", na.rm = TRUE, show.legend = FALSE, color = "#FF6F1B") +
                geom_line(aes(group = Ptid), alpha = 0.3, color = "#FF6F1B") +
                geom_boxplot(width = 0.25, alpha = 0.3, stat = "boxplot", outlier.shape = NA, show.legend = FALSE, color = "#FF6F1B") +
                # The lower and upper hinges correspond to the first and third quartiles (the 25th and 75th percentiles)
                # Whisker: Q3 + 1.5 IQR
                #scale_color_manual(name = "", values = c("#FF6F1B", "#0AB7C9"), guide = "none") + # guide = "none" in scale_..._...() to suppress legend
                geom_point(size = 3, alpha = 0.6, show.legend = TRUE, color = "#FF6F1B") +
                
                #geom_text(aes(label = ifelse(RespRate!="","Rate",""), x = 0.4, y = 5), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE) +
                geom_text(aes_string(x = x.var, label = "RespRate", y = 4.4), color = "black", size = panel.text.size, check_overlap = TRUE) +
                
                geom_hline(aes(yintercept = ifelse(RespRate!="",lbval,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
                geom_text(aes(label = ifelse(RespRate!="",lb,""), x = 0.4, y = lbval), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE, na.rm = TRUE) + 
                # only plot uloq for ID50
                geom_hline(aes(yintercept = ifelse(RespRate!="",lbval2,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
                geom_text(aes(label = ifelse(RespRate!="",lb2,""), x = 0.4, y = lbval2), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE, na.rm = TRUE) + 
                
                scale_x_discrete(labels = x.lb, drop=TRUE) +
                scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], 2), labels = scales::math_format(10^.x)) +
                labs(x = "Assay", y = unique(panel), title = paste(unique(panel), "longitudinal plots across timepoints"), color = "Category", shape = "Category") +
                plot_theme +
                guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
        
    return(p2)
}

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
            # write.csv(data.frame(x = x, y = y, strata = st, weight = wt), "input_columns.csv", row.names = FALSE)
            
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
            # saveRDS(corvec, file = "corvec.RDS")
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
    rr.x <- range(dat.tmp[, 1: (ncol(dat.tmp)-1)], na.rm = TRUE)
    rr.y <- range(dat.tmp[, 2: ncol(dat.tmp)], na.rm = TRUE)
    
    if (rr[1] == rr[2]) {
        rr <- c(rr[1] - 1, rr[2] + 1)
    }
    
    if (rr[2] - rr[1] < 2) {
        rr <- floor(rr[1]):ceiling(rr[2])
    }
    
    breaks <- floor(rr[1]):ceiling(rr[2])
    
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
        scale_x_continuous(limits = rr.x, breaks = breaks) + ylim(0, 1.25)
    
    for (j in 2:pairplots$nrow) {
        for (k in 1:(j - 1)) {
            pairplots[j, k] <- pairplots[j, k] +
                stat_smooth(
                    method = "loess", color = "red", se = FALSE,
                    lwd = loess_lwd
                ) +
                scale_x_continuous(
                    limits = rr.x, breaks = breaks,
                    labels = scales::math_format(10^.x)
                ) +
                scale_y_continuous(
                    limits = rr.y, breaks = breaks,
                    labels = scales::math_format(10^.x)
                )
        }
        pairplots[j, j] <- pairplots[j, j] +
            scale_x_continuous(
                limits = rr, # use maximum range of rr.x and rr.y here in order to show a complete histogram
                breaks = breaks,
                labels = scales::math_format(10^.x)
            ) + ylim(0, 1.25)
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


#' Pairplots of assay readout over time
#'
#' Produce the pairplots of assay readouts. The correlation is calculated by
#' the resampling-based strata adjusted Spearman rank correlation
#'
#' @param plot_dat: data frame: data for plotting.
#' @param assay: vector of strings: the assay name for plotting.
#' @param times: vector of strings: the time points for plotting.
#' @param strata: string: the column name in plot_dat that indicates the
#'  strata.
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
#'
#' @return pairplots: a ggplot object of the pairplot
covid_corr_pairplots_by_time <- function(plot_dat, ## data for plotting
                                         assay,
                                         times,
                                         strata,
                                         weight,
                                         plot_title,
                                         column_labels,
                                         height = 5.1,
                                         width = 5.05,
                                         units = "in",
                                         corr_size = 5,
                                         point_size = 0.5,
                                         loess_lwd = 1,
                                         plot_title_size = 10,
                                         column_label_size = 6.5,
                                         axis_label_size = 9,
                                         filename) {
    dat.tmp <- plot_dat[, paste0(times, assay)]
    rr <- range(dat.tmp, na.rm = TRUE)
    
    if (rr[1] == rr[2]) {
        rr <- c(rr[1] - 1, rr[2] + 1)
    }
    
    if (rr[2] - rr[1] < 2) {
        rr <- c(floor(rr[1]), ceiling(rr[2]))
    }
    
    breaks <- floor(rr[1]):ceiling(rr[2])
    
    if (rr[2] > ceiling(rr[1])) {
        breaks <- ceiling(rr[1]):floor(rr[2])
    } else {
        breaks <- floor(rr[1]):ceiling(rr[2]) ## breaks on the axis
    }
    
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
                     strata = subdat[, strata],
                     weight = subdat[, weight]
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
            panel.grid.minor = element_blank()
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
                    labels = label_math(10^.x)
                ) +
                scale_y_continuous(
                    limits = rr, breaks = breaks,
                    labels = label_math(10^.x)
                )
        }
        pairplots[j, j] <- pairplots[j, j] +
            scale_x_continuous(
                limits = rr, breaks = breaks,
                labels = label_math(10^.x)
            ) + ylim(0, 1.2)
    }
    
    ggsave(
        filename = filename, plot = pairplots, width = width, height = height,
        units = units
    )
}
