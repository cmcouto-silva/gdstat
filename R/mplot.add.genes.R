#' Title
#' 
#' Description
#' 
#' @usage mplot.add.genes(data, mplot, gene_col = "GENE", merge_col = "W_ID", by = "CHR", label = "geom_label_repel")
#' 
#' @param data description
#' @param mplot description
#' @param gene_col description
#' @param skip_empty_genes description
#' @param merge_col description
#' @param by description
#' @param label geom_text_repel or geom_label_repel
#' @param log_pvalue description
#' @param ... additional parameters for ggrepel functions
#' 
#' @return description
#' 
#' @import data.table ggplot2
#' @export

mplot.add.genes <- function(data, mplot, gene_col = "GENE", skip_empty_genes = T, merge_col = "W_ID",
                            by = "CHR", label = "geom_label_repel", log_pvalue = 2, ...) {
  
  # if(!any(grepl(col, colnames(data)))) {
  #   stop("Column name does not match with any column in data.")
  # }
  # 
  # if(!any(grepl(merge_col, colnames(data)))) {
  #   stop("Invalid column name passed on 'merge_col'.")
  # }
  
  if (!label %in% c("geom_text_repel", "geom_label_repel")){
    stop('Ony "geom_text_repel" and "geom_label_repel" can be specified in the label parameter.')
  }
  
  y.axis <- as.character(mplot$mapping$y[2])
  datagen <- merge(data, mplot$data[, .SD,, .SDcols = c(merge_col, "position")], by = merge_col, sort = F)
  
  if(skip_empty_genes) datagen <- datagen[GENE != ""]
  data.table::setcolorder(datagen, c("CHR", "POS", "SNP"))
  
  datagen <- datagen[datagen[, .I[which.max(get(y.axis))], by = by]$V1][LOG_PVALUE >= log_pvalue]
  
  # Plot genes
  if(label == "geom_text_repel") {
    plt <- mplot +
      ggrepel::geom_text_repel(data = datagen, mapping = aes(x = position, y = get(y.axis)),
                               label = datagen[, get(gene_col)], size = 2.5,
                               vjust = 1, nudge_x = 1)
  } else {
    plt <- mplot +
      ggrepel::geom_label_repel(data = datagen, mapping = aes(x = position, y = get(y.axis)),
                                label = datagen[, get(gene_col)], size = 3,
                                vjust = 1, nudge_x = 1)
  }
  
  plt$datagen <- datagen
  
  return(plt)
  
}
