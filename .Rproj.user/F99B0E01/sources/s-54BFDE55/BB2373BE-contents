#' Title
#' 
#' Description
#' 
#' @usage pbs_fst(geno.info, focal, close, outgroup, monomorphic_removal = T, method = "reynolds")
#' 
#' @param geno.info description
#' @param focal description
#' @param close description
#' @param outgroup description
#' @param monomorphic_removal description
#' @param method description
#' 
#' @return description
#' 
#' @import data.table
#' @export

pbs_fst <- function(geno.info, focal, close, outgroup, monomorphic_removal = T, method = "reynolds") {

  # Checking Arguments
  # if(class(geno.info) != "geno.info") stop("Data must be from pbs::geno.info class.")
  if(!method %in% c("reynolds")) stop('Only reynolds\' method is available.')
  if(!monomorphic_removal %in% c(T,F)) stop("Only logical values (TRUE or FALSE) are acceptable.")

  pops <- c(focal, close, outgroup)
  idx <- !pops %in% geno.info$frq[, CLST]

  if(any(idx)) stop("\n", paste0('Population "', pops[idx], '"', " did not find in CLST column.\n"))

  frq <- data.table(focal = geno.info$frq[CLST == focal, MAF],
                    close = geno.info$frq[CLST == close, MAF],
                    outgroup = geno.info$frq[CLST == outgroup, MAF]
  )

  # Removin monomorphic alleles in at least two populations
  if(monomorphic_removal == TRUE) {

    m0 <- frq[, lapply(.SD, function(.) . == 0L)][, I := .I]
    m0 <- m0[, M := sum(c(focal, close, outgroup)), by = I][, M >= 2L]

    m1 <- frq[, lapply(.SD, function(.) . == 1L)]
    m1 <- m1[, I := .I][, M := sum(c(focal, close, outgroup)), by = I][, M >= 2L]

    mono <- !(m0 | m1)
    cat(" ", sum(!mono), "monomorphic alleles were removed!\n")

    frq <- frq[mono]
    geno.info$bim <- geno.info$bim[mono]

  }

  # Applying FST calculation
  if (method == "reynolds") {

    fst.rwc <- function(p1, p2) {
      n <- ((p1 - p2)^2) + (((1 - p1) - (1 - p2))^2)
      d <- 2 * ( 1 - ((p1 * p2) + ((1 - p1) * (1 - p2))) )
      res <- n/d
      res[is.nan(res)] <- 0L
      return(res)
    }

    frq[, names(frq) := .(fst.rwc(focal, close), fst.rwc(focal, outgroup), fst.rwc(close, outgroup))]
    col_names <- c("focal.close","focal.outgroup", "close.outgroup")
    data.table::setnames(frq, names(frq), col_names)

  } else {
    stop ("At this moment, only Reynolds et al. (1983) formula has been implemented.")
  }

  # Returning data.table
  frq[, PBS := ( (-log(1L - focal.close)) + (-log(1L - focal.outgroup)) - (-log(1L - close.outgroup)) ) / 2L]
  frq[PBS < 0L, PBS := 0L] # set to 0 negative values

  return( cbind(geno.info$bim[, .(CHR, CM, POS, SNP)], frq[, .(PBS)]) )

}
