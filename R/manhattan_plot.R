#' Title
#' 
#' Description
#' 
#' @usage manhattan_plot(pbs_data, col_name, fig_name = NULL, facet_key = NULL)
#' 
#' @param pbs_data description
#' @param col_name description
#' @param fig_name description
#' @param facet_key description
#' @param point_size description
#' 
#' @return description
#' 
#' @import data.table
#' @export

manhattan_plot <- function(pbs_data, col_name, fig_name = NULL, facet_key = NULL, point_size = 0.5) {
  
  pbs_dt <- copy(pbs_data)
  pbs_dt[, position := .GRP, keyby = .(CHR, POS)]
  NCHR = pbs_dt[.N, CHR]
  
  pos_shift <- pbs_dt[, round(0.2 * .N / NCHR)]
  pbs_dt[, position := position + (CHR - 1) * pos_shift]
  # pbs_th <- quantile(pbs_dt[[col_name]], seq(0.995,0.999,0.001), na.rm = T)[c(1,5)]
  
  b_info <- pbs_dt[, .(breaks = 0.5 * (first(position) + last(position))), by = CHR]
  
  # Set colors
  plot_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  plot_colors <- rep_len(plot_colors, length(unique(pbs_dt[, CHR])))
  colors_vec <- as.factor(pbs_dt[, CHR])
  levels(colors_vec) <- plot_colors
  pbs_dt[, color := as.character(colors_vec)]
  
  # Plotting
  par(cex=2.5)
  gm <- ggplot2::ggplot(pbs_dt, aes_string(x = "position", y = col_name) ) + 
    geom_point(aes(color = color), size = point_size) +
    theme_bw(base_family = "arial", base_size = 12L) + 
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", size = 0.3),
          legend.position = "none") + 
    labs(x = "chromosome") +
    scale_x_continuous(breaks = unlist(b_info[, breaks]), 
                       label = unlist(b_info[, CHR])) #+
  # geom_hline(yintercept = pbs_th, lty = 2, col = gray(.1))
  
  if(!is.null(facet_key))
    gm <- gm +
    facet_wrap(as.formula(paste('~', facet_key)), 
               ncol = 1, scales = 'free')
  
  if(!is.null(fig_name))
    ggsave(fig_name)
  
  return(gm)
  
}
