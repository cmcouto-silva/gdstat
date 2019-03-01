#' Title
#' 
#' Description
#' 
#' @usage reynolds_fst(allele.freq, pop1, pop2)
#' 
#' @param allele.freq description
#' @param pop1 description
#' @param pop2 description
#' 
#' @return description
#' 
#' @export

reynolds_fst <- function(allele.freq, pop1, pop2) {
  
  if(missing(pop1) | missing(pop2))
    stop("Both populations must be specified!")
  
  allele.freq <- allele.freq$frq
  fst <- allele.freq[CLST == pop1, .(CHR, SNP, A1, A2)]
  
  p1 <- allele.freq[CLST == pop1, MAF] # rm [1:nsnp]
  p2 <- allele.freq[CLST == pop2, MAF] # rm [1:nsnp]
  
  n <- ((p1 - p2)^2) + (((1 - p1) - (1 - p2))^2)
  d <- 2 * ( 1 - ((p1 * p2) + ((1 - p1) * (1 - p2))) )
  
  d2 <- n/d
  d2[is.nan(d2)] <- 0L
  fst[, paste0("FST.", pop1, ".", pop2) := d2][]

}

# Example
# reynolds_fst(allele.freq = "plink.frq.strat", pop1 = "CEU", pop2 = "CHB")
