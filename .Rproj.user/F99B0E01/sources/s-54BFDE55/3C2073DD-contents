#' Return the index based in the number of SNPs for rolling statiscs
#'
#'\code{rolling_index_by_snp_number()} estimates the branch size for every SNP 
#' 
#' @import data.table
#' 
#' @param  focal_pop: data.table with the following columns: CHR, SNP, CM, POS,
#' NCHROBS, POP, VAR, AF.
#' @param close_pop: data.table
#' @param out_pop: data.table
#' 
#' @return data.table with the estimated branch size for every SNP

rolling_index_by_snp_number <-
  function(data_in, window_s = 3, step_s = 1){
  
    DT <- copy(data_in) 
    DT[, idx := .I]
    DT[, idd := .I]
    DT[, next_idd := idd + window_s - 1]
    setkey(DT, CHR, POS)

    if(step_s == 1){
      res <- DT[DT, .(idx = i.idx, seq = i.idx:idx),
                by = .EACHI,
                roll = Inf, 
                on = c(CHR = 'CHR', idd = 'next_idd')]
    }else{
      DT[, sid := seq(1:.N), by = CHR]
      DT[, ss := sid %% step_s ]
      res <- DT[DT, .(sid = i.sid, idx = i.idx, seq = i.idx:idx),
                by = .EACHI,
                roll = Inf, 
                on = c(CHR = 'CHR', idd = 'next_idd')]

      res <- res[sid %% step_s == 1]
    }
    return(list(res, DT))
  }
