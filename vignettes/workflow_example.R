# Set Working Directory
setwd("/home/cmcouto-silva/cmcouto.silva@usp.br/R/packages/gdstat")

# Load R scripts from package
devtools::document()
devtools::load_all()
# devtools::install_github("cmcouto-silva/gdstat")

# Load required libraries
# library(data.table)

#### Main Analysis ####
####################################################################################################!

# Common Arguments
plink <- "~/cmcouto.silva@usp.br/lab_files/datasets/SHGDP/unphased/chr2"
target_snps <- c("2:109513601", "11:61597212", "11:61580635")

focal <- "EastAsia"
close <- "WestEurasia"
outgroup <- "Africa"

# --- Calculate PBS --- #

# Calculate PBS with Fst estimates from Weir & Cockerham (1984)

pbs <- pbs_fst(plink = plink, focal = focal, close = close, outgroup = outgroup,
               monomorphic_removal = T, maf = 0.05, filter = "FID", method = "wc")

# Adding columns for genes
ref <- fread("~/cmcouto.silva@usp.br/lab_files/datasets/SHGDP/annot/shgdp_refGene.annot")
pbs <- merge(pbs, ref, by = c("CHR","POS"), sort = F)
pbs[, GENE := gsub("\\s*\\([^\\)]+\\)","", GENE)]

# # Exploratory Analysis
# diff(c(pbs[which(SNP == "2:109513601")-20L, POS],
#        pbs[which(SNP == "2:109513601")+20L, POS])
#      )
# 
# formatNumber <- function(n) {
#   n <- gsub('(?=(?:.{3})+$)', ".", "109513601", perl = TRUE)
#   n <- gsub("(^\\.)", "", n)
#   return(n)
# }

# --- Each SNP --- #

# Adding p-values
pbs[, P_VALUE := 1 - frankv(PBS, ties = 'random') / (.N + 1)]
# pbs[, P_VALUE := p.adjust(P_VALUE, "fdr", .N)]
pbs[, LOG_PVALUE := -log(P_VALUE, 10)]
pbs[LOG_PVALUE < 0, LOG_PVALUE := 0]
pbs[LOG_PVALUE >= 2]

# Getting peaks
pbs_peaks <- pbs[LOG_PVALUE >= 2]

# Plotting
mplot <- manhattan_plot(pbs_data = pbs, col_name = 'LOG_PVALUE')
mgplot <- mplot.add.genes(data = pbs_peaks, mplot = mplot, by = "CHR", merge_col = "SNP", label = "geom_label_repel")
mpplot <- mplot.add.points(DT = pbs, rolling_DT = NULL, mplot = mgplot, target_snp = target_snps, merge_col = "SNP", uniq = F)
mpplot

# --- SNP Mean --- #

# # Analyze PBS by sliding window mean
pbs_mean <- sliding_window(DT = pbs, window_size = 20, step = 5, mode = "fast", col = "PBS", keep.cols = T)

# Adding p-values
pbs_mean[, P_VALUE := 1 - frankv(PBS, ties = 'random') / (.N + 1)]
# pbs_mean[, P_VALUE := p.adjust(P_VALUE, "fdr", .N)]
pbs_mean[, LOG_PVALUE := -log(P_VALUE, 10)]
pbs_mean[LOG_PVALUE >= 2]

# Getting peaks
pbs_mean_peaks <- pbs_mean[LOG_PVALUE >= 2]
update_fields(DT = pbs, rolling_DT = pbs_mean_peaks, col = "LOG_PVALUE") # get biggest PBS and related-info by window
pbs_mean_peaks_nocons_windows <- exclude_consecutive_windows(pbs_mean_peaks) # exclude consecutive windows

pos <- pbs[, which(SNP %in% target_snps)]
pbs[pos:(pos+20)]

# RollingMean
mplot <- manhattan_plot(pbs_data = pbs_mean, col_name = 'LOG_PVALUE')
mgplot <- mplot.add.genes(data = pbs_mean_peaks_nocons_windows, mplot = mplot, by = "CHR", label = "geom_label_repel")
mpplot <- mplot.add.points(DT = pbs, rolling_DT = pbs_mean, mplot = mgplot, target_snp = target_snps, merge_col = "W_ID", uniq = F)
mpplot
