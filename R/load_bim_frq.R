#' Title
#' 
#' Description
#' 
#' @usage load_bim_frq(file_path, bim, frq)
#' 
#' @param file_path description
#' @param bim description
#' @param frq description
#' 
#' @return description
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
