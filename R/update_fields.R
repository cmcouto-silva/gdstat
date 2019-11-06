#' @export
update_fields <- function(DT, rolling_DT, col = "PBS", window_s = 20, step_s = 5, inplace = TRUE) {
  
  ref_idx <- refIdx(nrow(DT), window_s, step_s)
  
  for(i in rolling_DT[, W_ID]) {
    rolling_DT[W_ID == i, `:=`(
      setdiff(colnames(DT), c("PBS", "P_VALUE", "LOG_PVALUE")),
      DT[ref_idx[wid == i, seq(from, to)], .SD, .SDcols = setdiff(colnames(DT), c("PBS", "P_VALUE", "LOG_PVALUE"))][
        DT[ref_idx[wid == i, seq(from, to)], which.max(get(col))]]
    )]
  }
  return(rolling_DT)
}
