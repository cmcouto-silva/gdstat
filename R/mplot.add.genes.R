#' Title
#' 
#' Description
#' 
#' @usage mplot.add.genes(data, mplot, maxv_per_chr = T, all_peakv = T)
#' 
#' @param data description
#' @param mplot description
#' @param maxv_per_chr description
#' @param all_peakv description
#' 
#' @return description
#' 
#' @import data.table
#' @export

mplot.add.genes <- function(data, mplot, maxv_per_chr = T, all_peakv = T) {
  
  if(!any(grepl("GENE", colnames(data))))
    stop("GENE column not found.")
  
  if(maxv_per_chr) {
    
    # Getting max values per chromosome
    maxv_per_chr <- pbs_mean_peaks[pbs_mean_peaks[, .I[which.max(PBS)], by = CHR]$V1]
    setkey(mplot$data, CHR, CM, POS, SNP, PBS)
    setkey(maxv_per_chr, CHR, CM, POS, SNP, PBS)
    maxv_per_chr <- maxv_per_chr[mplot$data][!is.na(GENE)]
    
    # Plot genes
    mplot +
      ggrepel::geom_label_repel(data = maxv_per_chr, mapping = aes(x = position, y = PBS),
                                label = maxv_per_chr[, GENE], size = 3L, vjust = 1)
  }
  
  if(maxv_per_chr) {
    
    # Getting all peak values
    maxv_per_chr <- pbs_mean_peaks[pbs_mean_peaks[, .I[which.max(PBS)], by = CHR]$V1]
    setkey(mplot$data, CHR, CM, POS, SNP, PBS)
    setkey(maxv_per_chr, CHR, CM, POS, SNP, PBS)
    maxv_per_chr <- maxv_per_chr[mplot$data][!is.na(GENE)]
    
    # Plot genes
    mplot +
      ggrepel::geom_label_repel(data = maxv_per_chr, mapping = aes(x = position, y = PBS), 
                                label = maxv_per_chr[, GENE], size = 3L, vjust = 1)
  }
  
}
