#' Title
#' 
#' Description
#' 
#' @usage mplot.add.points(DT, rolling_DT = NULL, mplot, target_snp, merge_col = "W_ID", uniq = FALSE)
#' 
#' @param DT description
#' @param rolling_DT description
#' @param mplot description
#' @param merge_col description
#' @param uniq description
#' 
#' @return description
#' 
#' @import data.table ggplot2
#' @export

mplot.add.points <- function(DT, rolling_DT = NULL, mplot, target_snp, merge_col = "W_ID", uniq = FALSE) {
  
  if(is.null(rolling_DT)) {
    dataplot <- copy(DT)
    dataplot <- dataplot[SNP %in% target_snp]
  } else {
    snp.list <- list()
    rolling_DT = copy(rolling_DT)
    
    for(i in seq_along(target_snp)) {
      snp.list[[i]] <- search_target(DT, rolling_DT, target_snp[i])
      snp.list[[i]][, SNP := target_snp[i]][][]
    }
    
    rolling_DT <- rbindlist(snp.list)
    if(nrow(rolling_DT) == 0) warning("No SNP found in rolling_DT!")
    
    dataplot <- rolling_DT
  }
    
  y.axis <- as.character(mplot$mapping$y)[2]
  mpplot <- merge(dataplot, mplot$data[, .SD,, .SDcols = c(merge_col, "position")], by = merge_col, sort = F)
  
  if(uniq == T) {
    mpplot <- mpplot[mpplot[, .I[which.max(get(y.axis))], by = SNP]$V1] # unique = TRUE
  }
  
  mplot +
    ggplot2::geom_point(data = mpplot, mapping = aes(x = position, y = get(y.axis)))
  
}

