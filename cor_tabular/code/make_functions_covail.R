#' Generate response calls and fold-rise indicator for endpoints:
#' Responders at each pre-defined timepoint are defined as participants who
#' had baseline values below the LLOD with detectable ID80 neutralization titer
#' above the assay LLOD, or as participants with baseline values above
#' the LLOD with a 4-fold increase in neutralizing antibody titer.
#'
#' @param data The dataframe with the endpoint of interest at baseline and
#'  post-baseline.
#' @param bl The variable of endpoint value at baseline
#' @param post The variable of endpoint value at post-baseline
#' @param folds Folds to generate fold-rise indicator
#' @param cutoff LLOD or LLOQ
#' @param respoderFR Fold-rise used for response call positivity criteria
#'
#' @return Response calls and fold-rise indicator for \code{bl}: \emph{bl}Resp,
#'  \emph{bl}FR\emph{folds}
getResponder <- function(data,
                         # cutoff.name,
                         times=times, 
                         assays=assays, 
                         pos.cutoffs = pos.cutoffs) {
  
  # cutoff <- get(paste0("l", cutoff.name, "s"))
  for (i in times){
    for (j in assays){
      bl <- ifelse(substr(j, nchar(j)-1, nchar(j)) %in% c("la", "na", "sa"), paste0("B", substr(j, 1, nchar(j)-2)), paste0("B", j))
      post <- paste0(i, j)
      delta <- paste0("Delta", gsub("Day", "", i), "overB", j)
      cutoff <- pos.cutoffs[j]
      
      if(post %in% names(data)){
        data[, paste0(post, "Resp")] <- as.numeric(data[, post] > log10(cutoff))
        data[, delta] <- data[, post] - data[, bl]
      }
    }
  }
  return(data)
}



