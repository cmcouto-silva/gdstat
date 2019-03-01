#' Evaluete the sliding window meadian 
#'
#'\code{rolling_median_n()} evaluetes the rolling mean 
#' 
#' @import data.table
#' @export
#' 
#' @param data_in: data.table with the following columns: CHR, POS, col_name
#' @param col_name: name of the columns to be take the mean
#' @param widows_s: number of SNPs in the rolling window
#' @param step_s: number of SNPs in the rolling window
#' 
#' @return data.table with the rolling mean 

rolling_median_n <-
  function(data_in, col_name, window_s = 3, step_s = 2){

  aux <- 
    rolling_index_by_snp_number(data_in, 
                                window_s = window_s, 
                                step_s = step_s)

  res <- aux[[1]]
  DT <- aux[[2]]
  setnames(DT, col_name, 'aux')

  if(step_s == 1){
    DT[, MEDIAN := DT[res$seq, median(aux), by = res$idx]$V1]
    DT[, c('aux', 'idx', 'idd', 'next_idd') := NULL]
  }else{
    DT[sid %% step_s == 1, MEDIAN := DT[, .(DT[res$seq, median(aux), by = res$idx]$V1)]$V1]
    DT[, c('aux', 'idx', 'idd', 'next_idd', 'sid', 'ss') := NULL]
    DT <- DT[!is.na(MEDIAN)]
  }

  setnames(DT, 'MEDIAN', paste0(col_name, '_MEDIAN'))
  DT[, W_ID := .I]
  return(DT)
}
