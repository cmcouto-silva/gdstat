#' @export
update_fields <- function(DT, rolling_DT, col = "PBS", window_s = 20, step_s = 5, inplace = TRUE) {
  
  ref_idx <- refIdx(nrow(DT), window_s, step_s)
  
  for(i in rolling_DT[, W_ID]) {
    rolling_DT[W_ID == i, c("CHR", "SNP", "POS", "PBS") := DT[ref_idx[wid == i, seq(from, to)]][
      DT[ref_idx[wid == i, seq(from, to)], which.max(get(col))]][, .(CHR, SNP, POS, PBS)]
      ]
  }
}
