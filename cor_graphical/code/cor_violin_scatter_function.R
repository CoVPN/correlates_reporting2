
#' A ggplot object for violin box plot with or without lines
#' 
#' @param dat Dataframe with variables needed
#' @param dat.sample Random sample of the param dat for generating dots (showing all dots may be too much)
#' @param x X variable on x-axis
#' @param y Y variable on y-axis
#' @param colby a factor variable to specify box/dot/line/violin colors
#' @param shaby a factor variable to specify dot shapes
#' @param ylim Y-axis limits
#' @param ybreaks Y-axis breaks
#' @param ytitle X variable title
#' @param xtitle Y variable title
#' @param xlabel X variable label
#' @param toptitle Title for each page
#' @param type Type of figure: "noline" or "line"
#' @param facetby Faceting variables to form a matrix of panels
#' @param facetopt Faceting style: "wrap" or "grid"
#' @param group.num Number of case/non-case groups
#' @param col Color options for the colby param
#' @param shape Shape options for the shapeby param
#' @param col_lb Color breaks for the colby param
#' @param shp_lb Shape breaks for the shapeby param
#' @param prop.cex Font size for text within panels, n & response rate
#' @param ll.cex Font size for text within panels, eg: llod, pos.cut, uloq
#' @param rate.y.pos Y coordinate for showing response rate eg: "7.7"
#' @param pt.size point size
#' @param axis.text.x.cex font size for x axis text
#' @param axis.text.y.cex font size for y axis text
#' @param global.size global font size
#' @param n_rate variable for counts and response rate: "N_RespRate" or "N_RespRate_severe"
#' @return A ggplot object for violin + box plot with or without lines

violin_box_plot <- 
  function(dat, 
           dat.sample,
           x="time", 
           y="value", 
           colby="cohort_event", 
           shaby="cohort_event",
           ylim=c(1,7), 
           ybreaks=c(1,2,3,4,5,6),
           ytitle=NULL,
           xtitle="Time",
           xlabel=x_lb,
           toptitle=NULL,
           type="line",
           facetby=vars(cohort_event),
           facetopt="wrap",
           col=col_val,
           shape=shp_val,
           col_lb=cohort_event_lb,
           shp_lb=cohort_event_lb,
           prop.cex=5.4,
           group.num=3,
           ll.cex=prop.cex,
           rate.y.pos="7.7",
           n_rate,
           pt.size=5,
           axis.text.x.cex=25,
           axis.text.y.cex=25,
           global.size=11){
  
  p <- ggplot(data=dat, aes(x=.data[[x]], y=.data[[y]], color=.data[[colby]], shape=.data[[shaby]]))
  
  if (type=="line") {
    p <- p + geom_violin(scale="width", na.rm = TRUE)
      if (length(unique(dat.sample$time))!=1) p <- p + geom_line(data = dat.sample, aes(group = Ptid))
      # only draw line if there are multiple time points
      p <- p + geom_point(data = dat.sample, size = pt.size, show.legend = TRUE) +
      geom_boxplot(data = dat, width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)
  } else if (type=="noline") {
    p <- p + geom_violin(scale="width", na.rm = TRUE) +
      geom_jitter(data = dat.sample, width = 0.1, height = 0, size = pt.size, show.legend = TRUE) +
      geom_boxplot(data = dat, width=0.25, lwd=1.5, alpha = 0.3, outlier.shape=NA, show.legend = FALSE)
    }
  
  if (facetopt=="wrap") {p <- p + facet_wrap(facetby, ncol=group.num, drop=FALSE)
  } else if (facetopt=="grid") {p <- p + facet_grid(facetby, drop=FALSE)}
  
  p <- p + 
    geom_text(data = dat, aes(label=.data[[n_rate]], x=.data[[x]], y=as.numeric(rate.y.pos)), vjust = 1, color="black", size=prop.cex, check_overlap = TRUE) +
    geom_text(aes(label="n\nRate", x=0.4, y=rate.y.pos), vjust = 1, hjust = 0, color="black", size=prop.cex, check_overlap = TRUE) +
    geom_hline(aes(yintercept=lbval), linetype="dashed", color="gray", na.rm = TRUE) +
    geom_text(aes(label=lb, x=0.4, y=lbval), hjust = 0, color="black", size=ll.cex, check_overlap = TRUE, na.rm = TRUE) + 
    geom_hline(aes(yintercept=lbval2), linetype="dashed", color="gray", na.rm = TRUE) +
    geom_text(aes(label=lb2, x=0.4, y=lbval2), hjust = 0, color="black", size=ll.cex, check_overlap = TRUE, na.rm = TRUE) + 
    scale_x_discrete(labels=xlabel, drop=FALSE) +
    scale_y_continuous(limits=ylim, breaks=ybreaks, labels=math_format(10^.x)) +
    labs(x=xtitle, y=ytitle, title=toptitle, color="Category", shape="Category") +
    scale_color_manual(values=col, breaks=col_lb, drop=FALSE) +
    scale_shape_manual(values=shape, breaks=shp_lb, drop=FALSE) +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
          text = element_text(size=global.size),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=axis.text.x.cex),
          axis.text.y = element_text(size=axis.text.y.cex))

  return (p)
  }


#' A ggplot object for scatter plot
#' 
#' @param dat Dataframe with variables needed
#' @param x X variable on x-axis
#' @param y Y variable on y-axis
#' @param xtitle title for x-axis
#' @param ytitle title for y-axis
#' @param title main title
#' @param colby Variables to specify dot/line colors
#' @param shaby Variables to specify dot shapes
#' @param col Color options for the colby param
#' @param shape Shape options for the shapeby param
#' @param col_lb Color lbs for the colby param
#' @param shp_lb Shape lbs for the shapeby param
#' @param facetby facet variables
#' @param size point size
#' @param y.lim y-axis limits
#' @param y.breaks y-axis breaks
#' @param x.breaks x-axis breaks
#' @return A ggplot object for scatter plot

scatter_plot <- 
    function(dat, 
             x="Age", 
             y="value", 
             xtitle="Age (years)",
             ytitle=plots_ytitles[i],
             title=paste0(plots_titles[i],": ",timesls[[2]][d]),
             colby="cohort_event", 
             shaby="cohort_event",
             col=col_val,
             shape=shp_val,
             col_lb=gsub(" Cases", ifelse(case_set=="severe", " Severe Cases", " Cases"), cohort_event_lb),
             shp_lb=gsub(" Cases", ifelse(case_set=="severe", " Severe Cases", " Cases"), cohort_event_lb),
             facetby=as.formula("~Bserostatus+Trt"),
             size=with(dat, ifelse(cohort_event == "Non-Cases", 2.5, 4)),
             y.lim=c(floor(mins[plots[i]]), ceiling(maxs[plots[i]])), 
             y.breaks=seq(floor(mins[plots[i]]), ceiling(maxs[plots[i]])),
             x.breaks=seq(from=18, to=86, by=17)){
        
        p <- ggplot(dat, aes_string(x = x, y = y)) + 
            facet_wrap(facetby, nrow = 1) + 
            geom_point(alpha = 1, aes_string(color = colby, shape = shaby, size = size)) + 
            geom_smooth(aes_string(group = colby, color = colby), size=1.5, method = 'loess', se= F, span = 1.15) + 
            scale_y_continuous(limits=y.lim, breaks=y.breaks, labels=scales::math_format(10^.x)) +
            scale_x_continuous(breaks = x.breaks) +
            labs(title = title, x = xtitle, y = ytitle,
                 color="Category", shape="Category") +
            scale_color_manual(values = col, labels = col_lb, drop=F) +
            scale_shape_manual(values = shape, labels = shp_lb, drop=F) +
            guides(size = "none") +
            theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
                  panel.grid = element_blank(),
                  legend.title = element_text(size=22),
                  plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size=ifelse(c=="Vaccine", 27, 19)))
        
        return (p)
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
#' @param split.var group split variable in string, e.g., "panel", "assay_variant"
#' @param pointby a variable name by which different color and shape of point will be drawn, e.g, "cohort_col", "cohort_col2"
#' @param lgdbreaks breaks for point legend 
#' @param chtcols color panel for points
#' @param chtpchs shape panel for points
#' @param strip.text.y.size strip label size for y-axis, default is 25
#' @param axis.text.x.size x-axis label size, default is 9.5
#' @return A ggplot object list for longitudinal violin + box plot with lines

f_longitude_by_assay <- function(
  dat,
  x.var = "time_cohort",
  x.lb = c("BD1 Non-Cases","BD29 Non-Cases","BD1 Omicron Cases","BD29 Omicron Cases","DD1 Omicron Cases"),
  assays = assays,
  times = times,
  ylim = c(0,7.2),
  ybreaks = c(0,2,4,6),
  panel.text.size = 4,
  facet.x.var = vars(assay),
  facet.y.var = vars(Trt),
  split.var,
  pointby = "cohort_col",
  lgdbreaks = c("Omicron Cases", "Non-Cases", "Non-Responders"),
  chtcols = setNames(c("#FF6F1B", "#0AB7C9", "#8F8F8F"), c("Omicron Cases", "Non-Cases", "Non-Responders")),
  chtpchs = setNames(c(19, 19, 2), c("Omicron Cases", "Non-Cases", "Non-Responders")),
  strip.text.y.size = 25,
  axis.text.x.size = 9.5
) {
  
  plot_theme <- theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = axis.text.x.size),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 24, face="bold"),
          strip.text.x = element_text(size = 25), # facet label size
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
    mutate(cohort_col = ifelse(response==0 & !is.na(response), "Non-Responders", as.character(cohort_event)),
           time_cohort = factor(paste(time, cohort_event_adhoc),
                                levels = do.call(paste, expand.grid(c("Day 29", "Day 71", "Month 6"), c("Moderate Post-Peak Cases", "Severe Post-Peak Cases", "Non-Cases"))),
                                labels = do.call(paste, expand.grid(c("Day 29", "Day 71", "Month 6"), c("Moderate Post-Peak Cases", "Severe Post-Peak Cases", "Non-Cases"))))
    ) %>%
    ungroup() %>%
    group_split(.[[split.var]]) %>% # e.g., "panel" variable from assay_metadata
    purrr::map(function(d){
      ggplot(data = d, aes_string(x = x.var, y = "value", color = "cohort_event")) +
        facet_grid(rows = facet.y.var, col = facet.x.var) +
        
        geom_violin(scale = "width", na.rm = TRUE, show.legend = FALSE) +
        geom_line(aes(group = Ptid), alpha = 0.5) +
        geom_boxplot(width = 0.25, alpha = 0.3, stat = "boxplot", outlier.shape = NA, show.legend = FALSE) +
        # The lower and upper hinges correspond to the first and third quartiles (the 25th and 75th percentiles)
        # Whisker: Q3 + 1.5 IQR
        scale_color_manual(name = "", values = chtcols, guide = "none") + # guide = "none" in scale_..._...() to suppress legend
        # geoms below will use another color scale
        new_scale_color() +
        
        geom_point(aes(color = .data[[pointby]], shape = .data[[pointby]]), size = 3, alpha = 0.6, show.legend = TRUE) +
        scale_color_manual(name = "", values = chtcols, breaks = lgdbreaks, drop=FALSE) +
        scale_shape_manual(name = "", values = chtpchs, breaks = lgdbreaks, drop=FALSE) +
        
        geom_text(aes(label = ifelse(N_RespRate!="","Rate",""), x = 0.4, y = ylim[2] - 0.1), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE) +
        geom_text(aes_string(x = x.var, label = "N_RespRate", y = ylim[2] - 0.1), color = "black", size = panel.text.size, check_overlap = TRUE) +
        
        geom_hline(aes(yintercept = ifelse(N_RespRate!="",lbval,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
        geom_text(aes(label = ifelse(N_RespRate!="",lb,""), x = 0.4, y = lbval), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE, na.rm = TRUE) + 
        # only plot uloq for ID50
        geom_hline(aes(yintercept = ifelse(N_RespRate!="",lbval2,-99)), linetype = "dashed", color = "gray", na.rm = TRUE) +
        geom_text(aes(label = ifelse(N_RespRate!="",lb2,""), x = 0.4, y = lbval2), hjust = 0, color = "black", size = panel.text.size, check_overlap = TRUE, na.rm = TRUE) + 
        
        scale_x_discrete(labels = x.lb, drop=TRUE) +
        scale_y_continuous(limits = ylim, breaks = ybreaks, labels = scales::math_format(10^.x)) +
        labs(x = "Cohort", y = unique(d[[split.var]]), title = paste(unique(d[[split.var]]), "longitudinal plots across timepoints"), color = "Category", shape = "Category") +
        plot_theme +
        guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
    })
  return(p2)
}
