sliding_window_mean <- function(x, window_size = 20, slide = 5){
  idx1 <- seq(window_size, length(x), by = slide) + 1L
  idx2 <- seq(0L, by = slide, length.out = (length(idx1))) + 1L
  csum <- c(0L, cumsum(as.numeric(x)))
  
  return( (csum[idx1] - csum[idx2]) / window_size )
}

#' Title
#' 
#' Description
#' 
#' @usage sliding_window_mean_dt(dt, col_name, window_size = 20, slide = 5)
#' 
#' @param dt description
#' @param col_name description
#' @param window_size description
#' @param slide description
#' 
#' @return description
#' 
#' @export
sliding_window_mean_dt <- function(dt, col_name, window_size = 20, slide = 5){
  dt <- data.table::copy(dt)
  idx <- seq(window_size, nrow(dt), by = slide)
  idx1 <- idx + 1L
  idx2 <- seq(0L, by = slide, length.out = (length(idx1))) + 1L
  csum <- c(0L, cumsum(as.numeric(dt[[col_name]])))
  
  pbs_mean <- (csum[idx1] - csum[idx2]) / window_size
  dt <- dt[idx1]
  dt[, grep(col_name, colnames(dt), fixed = T) := pbs_mean][]

  return(dt)
}

