#' Evaluete the sliding window Moran statistics 
#'
#'\code{rolling_moran_n()} evaluetes the rolling mean 
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

rolling_moran_n <-
  function(data_in, col_name, window_s = 3, step_s = 2, qnorm_transform = TRUE){

  aux <- 
    rolling_index_by_snp_number(data_in, 
                                window_s = window_s, 
                                step_s = step_s)

  res <- aux[[1]]
  DT <- aux[[2]]
  setnames(DT, col_name, 'aux')

  moran <- 
    function(z){
      M <- z %*% t(z)
      n <- length(z)
      return((sum(M) - sum(diag(M))) / (n * n))
    }

  if(qnorm_transform){
    DT[, r := frankv(aux, ties = 'random') / (.N + 1)]
    DT[, z := qnorm(r)]
    DT[, r := NULL]
  }else{
    m <- DT[, mean(aux)]
    s <- DT[, sd(aux)]
    DT[, z := (aux - m) / s]
  }

  if(step_s == 1){
    DT[, MORAN := DT[res$seq, moran(z), by = res$idx]$V1]
    DT[, c('aux', 'idx', 'idd', 'next_idd', 'z') := NULL]
  }else{
    DT[sid %% step_s == 1, MORAN := DT[, .(DT[res$seq, moran(z), by = res$idx]$V1)]$V1]
    DT[, c('aux', 'idx', 'idd', 'next_idd', 'sid', 'ss', 'z') := NULL]
    DT <- DT[!is.na(MORAN)]
  }

  le <- log10(exp(1))
  DT[, p.value := - le * pnorm(MORAN, 0, 1, low = FALSE, log = TRUE)]

  setnames(DT, 'MORAN', paste0(col_name, '_MORAN'))
  DT[, W_ID := .I]
  return(DT)
}
