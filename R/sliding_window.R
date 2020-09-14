#' @export

sliding_window <- function(DT, window_size = 20, step = 5, mode = "fast", col = "PBS", keep.cols = T) {
  
  rolling_DT <- copy(DT)
  rolling_DT[, as.character(col) := frollmean(get(col), window_size, algo = mode, align = "left")]
  rolling_DT <- rolling_DT[seq(1, .N, by = step)][, W_ID := .I][!is.na(get(col))]
  if(isFALSE(keep.cols)) {
    rolling_DT <- rolling_DT[, .SD, .SDcols = c("CHR", "SNP", "POS", col, "W_ID")]
  }
  
  return(rolling_DT)
  
}
