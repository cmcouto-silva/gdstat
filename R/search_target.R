#' @export
search_target <- function(DT, rolling_DT, target_snp, window_s = 20, step_s = 5, inplace = FALSE) {
  
  ref_idx <- refIdx(nrow(DT), window_s, step_s)
  
  target_id <- DT[, which(SNP %in% target_snp)]
  target_wids <- ref_idx[target_id >= from & target_id <= to, wid]
  rolling_DT[W_ID %in% target_wids]
}
