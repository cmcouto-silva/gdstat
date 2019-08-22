#' @export
exclude_consecutive_windows <- function(rolling_DT, col = "LOG_PVALUE", inplace = F) {
  
  rolling_DT <- copy(rolling_DT)
  window_group <- c(diff(rolling_DT[, W_ID]) != 1)
  window_group <- rle(window_group)
  
  window_group_label <- c()
  for(i in seq_along(window_group$lengths))
    window_group_label <- append(window_group_label, rep(paste0("W", i), each = window_group$lengths[i]))
  
  W <- ifelse(rolling_DT[1, W_ID + 1] != rolling_DT[2, W_ID], "W0", "W1")
  window_group_label <- c(W, window_group_label)
  
  rolling_DT[, W_GROUP := window_group_label]
  rolling_DT[rolling_DT[, .I[which.max(get(col))], by = W_GROUP]$V1]
  
}
