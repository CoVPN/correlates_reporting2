#' A heatmap with annotation bars for biomarker correlation matrix.
#' 
#' @param dat_plot correlation matrix of markers of interest. We can get this by using cor().
#' @param main figure caption.
#' @param show_assay boolean. If TRUE, show assay color bar.
#' @param show_antigen boolean. If TRUE, show antigen color bar.
#' 
#' @return A heatmap of correlation.
#' 

cor_pheatmap <- function(
    dat_plot, 
    main, 
    size=10,
    show_assay = TRUE, 
    show_antigen = TRUE, 
    cluster_method = "euclidean",
    ...){
  
  #------------------------------
  # annotation for columns (markers)
  #------------------------------
  # Assay type
  Time <- sub(".*(B|Day28).*", "\\1", colnames(dat_plot))
  Time <- factor(Time,
                 levels = c("B", "Day28"),
                 labels = c("Baseline", "Peak"))
  
  # Assay type
  Antigen <- sub(".*[BDay28](FHA|FIM|PRN|PT|WCE).*", "\\1", colnames(dat_plot))
  Antigen <- factor(Antigen,
                  levels = c("FHA", "FIM", "PRN", "PT", "WCE"),
                  labels = c("FHA", "FIM", "PRN", "PT", "WCE"))
  
  # Assay type
  Assay <- groups <- case_when(
    str_detect(colnames(dat_plot), "Norm_Nasal_IgA$") ~ "Normalized Nasal IgA",
    str_detect(colnames(dat_plot), "Nasal_IgA$")      ~ "Nasal IgA",
    str_detect(colnames(dat_plot), "Serum_IgA$")      ~ "Serum IgA",
    str_detect(colnames(dat_plot), "Serum_IgG$")      ~ "Serum IgG")
  Assay <- factor(Assay,
                  levels = c("Normalized Nasal IgA", "Nasal IgA", "Serum IgA", "Serum IgG"),
                  labels = c("Normalized Nasal IgA", "Nasal IgA", "Serum IgA", "Serum IgG"))
  

  # annotations of columns
  annotations_col <- data.frame(
    Assay = Assay,
    Antigen = Antigen,
    Time = Time)
  rownames(annotations_col) <- colnames(dat_plot)
  colnames(annotations_col) <- c("Assay", "Antigen", "Time")
  
  
  
  #------------------------------
  # annotation colors
  #------------------------------
  # Define color palettes for heatmap and annotations
  # heatmap_col <- colorRampPalette(c("gray","yellow","green", "red", "black"))(50)
  heatmap_col <- colorRampPalette(c("white","royalblue4"))(50)
  Assay_col <- brewer.pal(length(levels(Assay)), "Set3")
  Antigen_col <- brewer.pal(length(levels(Antigen)), "Set1")
  
  # rename for seq control
  names(Assay_col) <- c("Normalized Nasal IgA", "Nasal IgA", "Serum IgA", "Serum IgG")
  names(Antigen_col) <- c("FHA", "FIM", "PRN", "PT", "WCE")
  

  # Example annotation colors
  annotation_colors <- list(
    Assay = Assay_col[levels(Assay)[levels(Assay) %in% unique(Assay)]],
    Antigen = Antigen_col[levels(Antigen)[levels(Antigen) %in% unique(Antigen)]],
    Time = c("Baseline" = "cornsilk", "Peak" = "gray25"))
  
  
  
  #--------------------
  # plotting
  #--------------------
  # color bar visibility
  if(!show_assay){
    annotations_col <- annotations_col[,!names(annotations_col) =="Assay"]
  }
  if(!show_antigen){
    annotations_col <- annotations_col[,!names(annotations_col) =="Antigen"]
  }
  
  # heatmap
  cor1 <- pheatmap(
    mat = dat_plot,
    scale = "none",
    # color palette for heatmap
    color = heatmap_col,
    # kmeans_k = 10,
    # display_numbers = TRUE,
    # whether if we want the trees
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # whether if we want to show row/colnames on the side
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_distance_rows = cluster_method,
    clustering_distance_cols = cluster_method,
    # annotation bars
    # annotation_row = annotations_col,
    annotation_col = annotations_col,
    # annotation colors
    annotation_colors = annotation_colors,
    # color of cell borders on heatmap, use NA if no border should be drawn
    border_color = NA,
    treeheight_row = 0,
    treeheight_col = 0,
    fontsize = size,
    fontsize_row = 15,
    fontsize_col = 15,
    # angle of the column labels: 90: (normal) counterclockwise 90
    angle_col = 90,
    # size of heatmap area
    # width = 10,
    # height = 10,
    # dendrogram ordering based on weights
    # clustering_callback = callback,
    # customize color keys
    # legend_breaks = c(-4, -3, -2, -1, 0, 1), # legend customisation
    # legend_labels = c("1e-4", "1e-3", "1e-2", "1e-1", "0", "1"), # legend customisation
    main = main,
    silent = TRUE
  )
  return(cor1)
}
