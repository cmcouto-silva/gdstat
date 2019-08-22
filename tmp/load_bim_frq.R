#' @title Load Plink files (.bim and .strat.frq)
#' 
#' @description Load SNP and allele frequency Plink data (.bim and .strat.frq, respectivelly), and save them as a list 
#' which is required for running further functinons.
#' 
#' @usage load_bim_frq(file_path, bim, frq)
#' 
#' @param file_path File path to file (without extensions).
#' @param bim Scalar character. File path to .bim file. 
#' @param frq Scalar character. File path to .strat.frq file.
#' 
#' @return list with SNP abd allele frequecy data.
#' 
#' @export

load_bim_frq <- function(file_path, bim, frq){
  
  if(!missing(bim)) {
    bim <- data.table::fread(bim, col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  } else {
    bim <- data.table::fread(paste0(file_path, '.bim'), col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  }
  
  if(!missing(frq)) {
    frq <- data.table::fread(frq)
  } else {
    frq <- data.table::fread(paste0(file_path, '.frq.strat'))
  }
    
  list(bim = bim, frq = frq)
  
}
