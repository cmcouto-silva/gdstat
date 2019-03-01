#' Title
#' 
#' Description
#' 
#' @usage get_peaks(data, quantile = 0.999, snp.annot = TRUE, rm.loc = TRUE)
#' 
#' @param data description
#' @param quantile description
#' @param snp.annot description
#' @param rm.loc description
#' 
#' @return description
#' 
#' @import data.table
#' @export

get_peaks <- function(data, quantile = 0.999, snp.annot = TRUE, rm.loc = TRUE) {

  # Managing arguments
  data <- copy(data)
  qtl <- quantile(x = data[, PBS], quantile)

  # Getting desired peaks
  data <- data[data[, .I[PBS >= qtl], by = CHR]$V1]

  # SNP Annotation
  if(snp.annot) {
    data[, GENE := snp.annot(SNP)][]
    setcolorder(data, c("CHR","CM","POS","SNP","GENE","PBS"))
    
    if(rm.loc) {
      data[, GENE := sapply(GENE, function(.) {
        genes <- unlist(strsplit(., ","))
        genes <- genes[!grepl("LOC", genes)]
        paste(genes, collapse = ",")
      })]
    }
    
    data[GENE == "", GENE := NA][]
  }

}
