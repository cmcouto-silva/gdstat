#' @export

sliding_window <- function(DT, window_size = 20, step = 5, mode = "fast", col = "PBS", FUN = "mean", keep.cols = F) {
  
  
  rolling_DT <- copy(DT)
  
  if (mode == "fast") {
    rolling_DT[, as.character(col) := frollmean(PBS, window_size, align = "left")]
    rolling_DT <- rolling_DT[seq(1, .N, by = step)][,W_ID := .I][!is.na(get(col))]
    if(isFALSE(keep.cols)) {
      rolling_DT <- rolling_DT[, .SD, .SDcols = c("CHR", "SNP", "POS", col, "W_ID")]
    }
  }
  
  if(mode == "exact") {
    
    ref_idx <- refIdx(nrow(rolling_DT), window_size, step)
    
    rolling_DT <- rbindlist(lapply(1:nrow(ref_idx), function(i) {
      window <- ref_idx[i, seq(from, to)]
      target_line <- rolling_DT[window, which.max(get(col))]
      rolling_DT[window, c(rolling_DT[target_line, .SD, .SDcols = setdiff(colnames(rolling_DT), col)], 
                           lapply(.SD, FUN)), .SDcols = col]
      }))
  }
  
  return(rolling_DT)
  
}

