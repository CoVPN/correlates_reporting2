
#' A ggplot object for violin box plot with or without lines
#' 
#' @param dat Dataframe with variables needed
#' @param dat.sample Random sample of the param dat for generating dots (showing all dots may be too much)
#' @param x X variable on x-axis
#' @param y Y variable on y-axis
#' @param colby Variables to specify box/dot/line/violin colors
#' @param shaby Variables to specify dot shapes
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
#' @param col_lb Color lbs for the colby param
#' @param shp_lb Shape lbs for the shapeby param
#' @param prop.cex Font size for text within panels, n & response rate
#' @param ll.cex Font size for text within panels, eg: llod, pos.cut, uloq
#' @param rate.y.pos Y coordinate for showing response rate eg: "7.7"
#' @param pt.size point size
#' @param axis.text.x.cex font size for x axis text
#' @param axis.text.y.cex font size for y axis text
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
           axis.text.y.cex=25){
  
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
    scale_color_manual(values=col, labels=col_lb, drop=FALSE) +
    scale_shape_manual(values=shape, labels=shp_lb, drop=FALSE) +
    theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "in"),
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