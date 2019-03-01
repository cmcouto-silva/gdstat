#' Estimate the FST using the Wier-Cockerhan formula  
#'
#'\code{wc_fst()} estimates the Fst for every variant 
#' 
#' @import data.table
#' @export
#' 
#' @param  pop_list: list of where each element is a data.table with the
#' following columns: CHR, SNP, CM, POS, NCHROBS, POP, VAR, AF
#' 
#' @return data.table with Fst estimated for every SNP 



wc_fst <-
    function(pop_list){

    r <- length(pop_list)
    pop_dt <- rbindlist(pop_list)

    setkeyv(pop_dt, c("CHR", "CM", "POS", "VAR", "SNP"))

    pop_dt <- pop_dt[pop_dt[, .N , by = .(CHR, CM, POS, VAR, SNP)][N == r]]
    pop_dt[, N := NULL]
    pop_dt[, G := .GRP, by = .(CHR, CM, POS, VAR, SNP)]

    # Assuming biallelic data select one variant to estimate the Fst 
    pop_dt <- pop_dt[G %% 2 != 0]
    pop_dt[, `:=`(G = NULL, VAR = NULL)]

    setkeyv(pop_dt, c("CHR", "CM", "POS", "SNP"))

    aux1 <- pop_dt[pop_dt[, .(nbar = mean(NCHROBS), 
                              nc = (sum(NCHROBS) - sum(NCHROBS ^ 2) / sum(NCHROBS)) / (r - 1), 
                              pbar = sum(NCHROBS * AF) / sum(NCHROBS)),
                          by = .(CHR, CM, POS, SNP)]
                  ]
    aux1[, G := .GRP, POP]

    aux2 <- 
      aux1[, .(ssqr = sum(NCHROBS * (pbar - AF) ^ 2) / ((r - 1) * nbar),
              pqbar = pbar * (1 - pbar), 
              nbar, 
              nc, 
              G
               ),
          by = .(CHR, CM, POS, SNP)
          ][G == 1]

    aux2[, G := NULL]

    fst_dt <- 
      aux2[,`:=`(T1 = ssqr - (pqbar - (r - 1) * ssqr / r) / (nbar - 1), 
                T2 = (nc - 1) * pqbar / (nbar - 1) + 
                     (1 + (r - 1) * (nbar - nc)/ (nbar - 1)) * ssqr / r)
          ][,
           `:=`(FST = T1 / T2)]

    fst_dt[, c("ssqr", "pqbar", "nbar", "nc") := NULL]

    fst_dt[ (FST < 0) | is.na(FST) , FST := 0] 

    return(fst_dt[])
    }

