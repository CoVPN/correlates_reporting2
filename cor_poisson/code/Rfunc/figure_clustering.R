#' A heatmap with annotation bars for biomarker clustering.
#' 
#' @param dat_plot input data.frame of markers of interest.
#' @param main figure caption.
#' @param size font size.
#' @param show_assay boolean. If TRUE, show assay color bar.
#' @param show_antigen boolean. If TRUE, show antigen color bar.
#' @cluster_method the way to measure similarity, default is the euclidean distance.
#' 
#' @return A heatmap of correlation.
#' 

# dat_plot = dat_all %>% select(starts_with("B"))
# main = "Clustering Heatmap - Baseline Markers"
# size=6
# show_assay = TRUE
# show_antigen = TRUE
# cluster_method = "euclidean"
fig_clustering <- function(
    dat_plot, main, 
    size = 12,
    show_assay = TRUE, 
    show_antigen = TRUE,
    cluster_method = "euclidean"){
  
  #------------------------------
  # annotation for rows (participants)
  #------------------------------
  # nasal treatment
  Trt_nasal <- dat_plot
  Trt_nasal <- dat_plot$Trt_Nasal
  Trt_nasal <- factor(Trt_nasal, 
                  levels = c(1, 0),
                  labels = c("Yes", "No"))
  
  # tdap treatment
  Trt_tdap <- dat_plot$Trt_tdap
  Trt_tdap <- factor(Trt_tdap, 
                  levels = c(1, 0),
                  labels = c("Yes", "No"))
  
  # Gender
  Gender <- dat_plot$Gender
  Gender <- factor(Gender,
                   levels = c("Male", "Female"),
                   labels = c("Male", "Female"))
  
  # Trial
  Trial <- dat_plot$Trial
  Trial <- factor(Trial,
                  levels = c("201", "202"),
                  labels = c("IB-201P", "IB-202P"))
  
  # annotations of columns
  annotations_row <- data.frame(Trt_tdap, Trt_nasal, Gender, Trial)
  rownames(annotations_row) <- rownames(dat_plot)
  colnames(annotations_row) <- c("Trt-Tdap", "Trt-Nasal", "Gender", "Trial")
  
  
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
  heatmap_col <- colorRampPalette(c("white", "royalblue4"))(50)
  Assay_col <- brewer.pal(length(levels(Assay)), "Set3")
  Antigen_col <- brewer.pal(length(levels(Antigen)), "Set1")
  
  # rename for seq control
  names(Assay_col) <- c("Normalized Nasal IgA", "Nasal IgA", "Serum IgA", "Serum IgG")
  names(Antigen_col) <- c("FHA", "FIM", "PRN", "PT", "WCE")
  
  
  # Example annotation colors
  annotation_colors <- list(
    Assay = Assay_col[levels(Assay)[levels(Assay) %in% unique(Assay)]],
    Antigen = Antigen_col[levels(Antigen)[levels(Antigen) %in% unique(Antigen)]],
    Time = c("Baseline" = "cornsilk", "Peak" = "gray25"),
    `Trt-Tdap` = c("Yes" = "chartreuse3", "No" = "brown2"),
    `Trt-Nasal` = c("Yes" = "chartreuse3", "No" = "brown2"),
    Gender = c("Male" = "lightskyblue", "Female" = "palevioletred1"),
    Trial = c("IB-201P" = "gold", "IB-202P" = "darkslategrey")
    )
  
  

  p1 <- pheatmap(
    mat = dat_plot[, !(colnames(dat_plot) %in% c("Trial", "Gender", "Trt_Nasal", "Trt_tdap"))],
    scale = "none",
    # color palette for heatmap
    color = heatmap_col,
    # kmeans_k = 10,
    # whether if we want the trees
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    # whether if we want to show row/colnames on the side
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_distance_rows = cluster_method,
    clustering_distance_cols = cluster_method,
    # annotation bars
    annotation_row = annotations_row,
    annotation_col = annotations_col,
    # annotation colors
    annotation_colors = annotation_colors,
    # color of cell borders on heatmap, use NA if no border should be drawn
    border_color = NA,
    treeheight_row = 50,
    treeheight_col = 50,
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
    # legend_labels = c("1e-4", "1e-3", "1e-2", "1e-1", "0", "10"), # legend customisation
    main = main,
    silent = TRUE
  )
  return(p1)
}

# 
# # antigen label
# Antigen_A <- sub(".*A([0-9]+).*", "\\1", colnames(dat_plot))
# Antigen <- c(paste(Type, Antigen_A, sep = "_"))
# Antigen <- factor(Antigen, 
#                   levels = c(paste0("T4_", sprintf("%02d", 1:13)),
#                              paste0("T8_", sprintf("%02d", 1:13)),
#                              paste0("IGG_", sprintf("%02d", 14:23)),
#                              paste0("BA_", sprintf("%02d", 45))),
#                   labels = c("CD4 E or M BA.1-2-4-5",
#                              "CD4 N BA.4-5",
#                              "CD4 S1 BA.4-5",
#                              "CD4 S2 BA.4-5",
#                              "CD4 S1 CON",
#                              "CD4 S2 CON",
#                              "CD4 E-M Index",
#                              "CD4 N Index",
#                              "seb control",
#                              "CD4 Any BA.4-5",
#                              "CD4 S BA.4",
#                              "CD4 Any Index",
#                              "CD4 S Index", # T4: 1-13
#                              "CD8 E or M BA.1-2-4-5",
#                              "CD8 N BA.4-5",
#                              "CD8 S1 BA.4-5",
#                              "CD8 S2 BA.4-5",
#                              "CD8 S1 CON",
#                              "CD8 S2 CON",
#                              "CD8 E-M Index",
#                              "CD8 N Index",
#                              "seb control",
#                              "CD8 Any BA.4-5",
#                              "CD8 S BA.4",
#                              "CD8 Any Index",
#                              "CD8 S Index", # T8: 1-13
#                              "IgG S BA.2.12.1",
#                              "IgG N Index", 
#                              "IgG S Index",
#                              "IgG S Alpha",
#                              "IgG S Beta",
#                              "IgG S Delta",
#                              "IgG S BA.1.1.529; BA.1",
#                              "IgG S BA.2",
#                              "IgG S BA.2.75",
#                              "IgG S BA.4/5", # IGG: 14-23
#                              "NAb BA.4/5")) 