#' @export
exclude_consecutive_windows <- function(rolling_DT, col = "LOG_PVALUE", inplace = F) {
  
  rolling_DT <- copy(rolling_DT)
  rolling_DT[, c("W_GROUP", "I") := .(character(), .I)]
  
  window_groups <- split(rolling_DT$W_ID, cumsum(c(1, diff(rolling_DT$W_ID) != 1)))
  
  for(i in seq_along(window_groups)) {
    set(rolling_DT, i = rolling_DT[W_ID %in% window_groups[[i]], I], j = "W_GROUP", value = paste0("WD_", i))
  }
  
  rolling_DT[, I := NULL]
  rolling_DT[rolling_DT[, .I[which.max(get(col))], by = W_GROUP]$V1]
  
}
