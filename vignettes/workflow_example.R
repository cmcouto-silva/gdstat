# Set Working Directory
setwd("/home/cmcouto-silva/cmcouto.silva@usp.br/R/packages/gdstat")

# Load R scripts from package
devtools::document()
devtools::load_all()

# Load required libraries
# library(data.table)

#### Main Analysis ####
####################################################################################################!

# Common Arguments
plink <- "inst/extdata/chr2_afr_eas_eur_maf"

focal <- "EUR"
close <- "EAS"
outgroup <- "AFR"

target_snp <- 136608646

# --- Calculate PBS --- #

# Calculate PBS with Fst estimates from Reynolds, Weir & Cockerham (1983)
pbs_rwc <- pbs_fst(plink = plink, focal = focal, close = close, outgroup = outgroup,
                   monomorphic_removal = T, filter = "FID", method = "rwc")

pbs_rwc <- pbs_rwc[!(is.na(focal.close) | is.na(focal.outgroup) | is.na(close.outgroup))]

# Calculate PBS with Fst estimates from Weir & Cockerham (1984)
pbs_wc <- pbs_fst(plink = plink, focal = focal, close = close, outgroup = outgroup,
                  monomorphic_removal = F, filter = "FID", method = "wc")

pbs_wc <- pbs_wc[!(is.na(focal.close.fst) | is.na(focal.outgroup.fst) | is.na(close.outgroup.fst))]

# --- PBS obtained by Reynolds, Weir & Cockerham (1983)'s Fst --- #

# Calculate p-value
pbs_rwc[, P_VALUE := 1 - frankv(PBS, ties = 'random') / (.N + 1)]
pbs_rwc[, LOG_PVALUE := -log(P_VALUE)]

# Calculate distribution of quantils
PBS_wc_distrQT <- quantile(pbs_wc[, PBS], c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99, 0.999))

# Verificar se o SNP-alvo é outlier
pbs_rwc[POS == target_snp]
pbs_wc[POS == target_snp]

fst_EurEas_distrQT <- quantile(pbs_wc[, focal.close.fst],  c(.01, .05, .1, .25, .50,  .75, .90, .95, .99))
fst_EurAfr_distrQT <- quantile(pbs_wc[, focal.outgroup.fst],  c(.01, .05, .1, .25, .50,  .75, .90, .95, .99))
fst_EasAfr_distrQT <- quantile(pbs_wc[, close.outgroup.fst],  c(.01, .05, .1, .25, .50,  .75, .90, .95, .99))

SNPfrom_BP <- target_snp - 10000
SNPto_BP <- target_snp + 10000

SNPfrom_id <- max(which(pbs_wc[, POS] <= SNPfrom_BP))
SNPto_id <- min(which(pbs_wc[, POS] >= SNPto_BP))

pbs_wc <- pbs_wc[SNPfrom_id:SNPto_id]

pbs_wc[PBS < 0]
pbs_wc[is.na(PBS)]

# EUR vs AFR
plot( x = pbs_wc[, POS], y = pbs_wc[, focal.outgroup.fst],
  xlab = "position", ylab = "Fst", main = "EUR vs AFR",
  ylim = c(0,1), pch = 20, cex = 0.5 )

points(x = pbs_wc[POS == target_snp, POS], pbs_wc[POS == target_snp, focal.outgroup.fst], col = "blue")
abline(h = fst_EurAfr_distrQT[["95%"]])

# Plot PBS
plot(x = pbs_wc[, POS], y = pbs_wc[, PBS], pch = 20, cex = 0.5)
abline(h = c(PBS_wc_distrQT["95%"], PBS_wc_distrQT["99%"]), lty = 2)
# abline(h = PBS_rwc_distrQT["99.9%"], lty = 2)
points(x = pbs_wc[POS == target_snp, POS], pbs_wc[POS == target_snp, PBS], col = "blue")

# Plot Fst

# EUR vs EAS
plot(x = pbs_rwc[range_snps, POS], y = pbs_rwc[range_snps, focal.close], pch = 20, cex = 0.5,
     main = "EUR vs EAS", xlab = "position", ylab = "Fst")
abline(h = quantile(x = pbs_rwc[, focal.close], c(0.95, 0.99)), lty = 2)
points(x = pbs_rwc[POS == target_snp, POS], pbs_rwc[POS == target_snp, focal.close], col = "blue")

# EUR vs AFR
plot(x = pbs_rwc[range_snps, POS], y = pbs_rwc[range_snps, focal.outgroup], pch = 20, cex = 0.5,
     main = "EUR vs AFR")
abline(h = quantile(x = pbs_rwc[range_snps, focal.outgroup], c(0.95, 0.99)), lty = 2)
points(x = pbs_rwc[POS == target_snp, POS], pbs_rwc[POS == target_snp, focal.outgroup], col = "blue")

# EAS vs AFR
plot(x = pbs_rwc[range_snps, POS], y = pbs_rwc[range_snps, close.outgroup], pch = 20, cex = 0.5,
     main = "EAS vs AFR")
abline(h = quantile(x = pbs_rwc[, close.outgroup], c(0.95, 0.99)), lty = 2)
points(x = pbs_rwc[POS == target_snp, POS], pbs_rwc[POS == target_snp, close.outgroup], col = "blue")





# Plot individual PBS values
png(paste0(gene, "/", method, "/", "pbs_", m[i,1], "_", m[i,2], "_", m[i,3], ".png"), width = 800, res = 100)
plot(pbs[, POS], pbs[, PBS], xlab = "POS", ylab = "PBS", main = paste(m[i, 1], "vs", m[i, 2], "vs", m[i, 3]))
points(pbs[POS %in% target_pos, POS], pbs[POS %in% target_pos, PBS], pch = 16, col = "#540907")
points(pbs[POS %in% target_pos_specific, POS], pbs[POS %in% target_pos_specific, PBS], pch = 16, col = "red")
abline(h = PBS_distrQT["99.9%"], col = 'red')
dev.off()



# PBS for individuals SNPs
pbs <- pbs_fst(geno.info = bim.frq, focal = "EUR", close = "EAS", outgroup = "AFR", monomorphic_removal = F)
# manhattan_plot(pbs_data = pbs, col_name = 'PBS', fig_name = 'pbs_whole.png')

PBS_distrQT <- quantile(pbs[, PBS], c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
PBS_distrQT # estimar os quantils

LCT <- 136608646
pbs[PBS > PBS_distrQT["99%"]]

lct <- ..

p_value_lct <- pbs[PBS >= pbs[POS == LCT, PBS], .N] / nrow(pbs)




p_values <- sapply(pbs[, PBS], function(p){
  pbs[PBS >= p, .N] / nrow(pbs)
})




lct_pos <- pbs[, which(POS == LCT)]





reynolds_fst(bim.frq, "AFR", "EAS")[lct_pos]
reynolds_fst(bim.frq, "AFR", "EUR")[lct_pos]
reynolds_fst(bim.frq, "EAS", "EUR")[lct_pos]

pbs

# PBS for mean of SNPs
pbs_mean <- sliding_window_mean_dt(pbs, 'PBS')
pbs_mean_peaks <- get_peaks(pbs_mean, quantile = 0.999, snp.annot = F, rm.loc = F)

pbs_mean[, P_VALUE := 1 - frankv(PBS, ties = 'random') / (.N + 1)]


lct_pos <- 136545415:136594750
136608646
snpinfo <- gt::annot_snp_genes(pbs[POS %in% lct_pos, SNP])


pbs_mean_peaks

window_size = 20
slide = 5




seq(window_size, 948, by = slide)

seq(1:948, by = 20, each = 5)



seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 



seqle(wid)
improved.rle(x = wid)


# Plotting
mplot <- manhattan_plot(pbs_data = pbs_mean, col_name = 'PBS')

mplot_peaks <- function(data, mplot, maxv_per_chr = T, all_peakv = T)


maxv_per_chr <- pbs_mean_peaks[pbs_mean_peaks[, .I[which.max(PBS)], by = CHR]$V1]
setkey(mplot$data, CHR, CM, POS, SNP, PBS)
setkey(maxv_per_chr, CHR, CM, POS, SNP, PBS)
maxv_per_chr <- maxv_per_chr[mplot$data][!is.na(GENE)]

all_peakv <- copy(pbs_mean_peaks)
setkey(all_peakv, CHR, CM, POS, SNP, PBS)
setkey(mplot$data, CHR, CM, POS, SNP, PBS)
all_peakv <- all_peakv[mplot$data][!is.na(GENE)]
all_peakv <- all_peakv[all_peakv[, .I[which.max(PBS)], by = GENE]$V1]

# Plot only max peak values by chromosome
mplot +
  ggrepel::geom_label_repel(data = maxv_per_chr, mapping = aes(x = position, y = PBS),
                            label = maxv_per_chr[, GENE], size = 3L, vjust = 1)

# Plot all peak values
mplot +
  ggrepel::geom_text_repel(data = all_peakv, mapping = aes(x = position, y = PBS),
                            label = all_peakv[, GENE], size = 2.5,
                           vjust = 0.5, nudge_x = -0.35)


uniq_genes <- all_peakv[, unique(GENE)]
uniq_genes <- uniq_genes[uniq_genes != ""]

pbs_mean <- rolling_mean_n(pbs, col_name = 'PBS', window_s = 20, step_s = 5)
pbs_mean[, P_RANK := 1 - frankv(PBS_MEAN) / (.N + 1)]
manhattan_plot(pbs_data = pbs_mean, col_name = 'PBS_MEAN', fig_name = 'pbs_mean.png')

pbs_p99 <- pbs_mean[PBS_MEAN >= quantile(PBS_MEAN, 0.999)]
pbs_p99[, GENE := LD::snp.annot(SNP)]
pbs_p99[CHR == 22]

# PBS for median of SNPs
pbs_median <- rolling_median_n(pbs, col_name = 'PBS', window_s = 20, step_s = 5)
manhattan_plot(pbs_data = pbs_median, col_name = 'PBS_MEDIAN', fig_name = 'pbs_median.png')

#### Analysing P-values ####
####################################################################################################!

# Individual SNPs
pbs[, P_RANK := 1 - frankv(PBS, ties = 'random') / (.N + 1)]
pbs[, LOG_PVAL := -log(P_RANK)]

# PBS for mean of SNPs
pbs_mean[, P_RANK := 1 - frankv(PBS_MEAN) / (.N + 1)]
pbs_mean[, LOG_PVAL := - log(P_RANK)]

# PBS for median of SNPs
pbs_median[, P_RANK := 1 - frankv(PBS_MEDIAN, ties = 'random') / (.N + 1)]
pbs_median[, LOG_PVAL := -log(P_RANK)]

# Getting Peaks
peaks_mean <- get_peaks(pbs, pbs_mean, stat_col = 'PBS',
                        score_col = 'P_RANK', score_th = 0.01,
                        window_s = 20, greater = FALSE)
peaks_median <-
  get_peaks(pbs, pbs_median, stat_col = 'PBS',
            score_col = 'P_RANK', score_th = 0.01,
            window_s = 20, greater = FALSE)


pbs_median[peaks_median[CHR == 2], on = .(CHR, CM, W_ID)]
pbs_median[, summary(LOG_PVAL)]


#### Extra section: comparing fst values ####
####################################################################################################!

# from Kelly
kelly_res <- fread("FST_AMZ2_MES_EAS.txt")

# herein
bim.frq <- load_bim_frq2(file_path = "nam_phased_adj")

res <- data.table (
  reynolds_fst(bim.frq, pop1 = "MES", pop2 = "BRZ")[, c(2,5)],
  reynolds_fst(bim.frq, pop1 = "EAS", pop2 = "BRZ")[, 5L],
  reynolds_fst(bim.frq, pop1 = "EAS", pop2 = "MES")[, 5L]
)

# MES NAM
merge(x = res[, .(SNP, FST.MES.BRZ)], y = kelly_res[, .(rs, FST.MES.NAM)],
      by.x = "SNP", by.y = "rs", sort = F)[, round(.SD, 4L), .SDcols = 2:3
                                           ]

# EAS NAM
merge(x = res[, .(SNP, FST.EAS.BRZ)], y = kelly_res[, .(rs, FST.EAS.NAM)],
      by.x = "SNP", by.y = "rs", sort = F)[, round(.SD, 4L), .SDcols = 2:3
                                           ]

# EAS MES
merge(x = res[, .(SNP, FST.EAS.MES)], y = kelly_res[, .(rs, FST.EAS.MES)],
      by.x = "SNP", by.y = "rs", sort = F)[, round(.SD, 4L), .SDcols = 2:3

                                                                                ]
#### Testing FST/PBS Mean ####
####################################################################################################!

# Loading scripts
devtools::load_all()

# Load .bim / .frq.strat files
bim.frq <- load_bim_frq2(file_path = "nam_phased_adj")

### PBS Mean
pbs_mean <- pbs_fst2(frq = bim.frq, focal = "BRZ", close = "MES", outgroup = "EAS")
pbs_mean <- sliding_window_mean_dt(pbs_mean, "PBS")
manhattan_plot(pbs_data = pbs_mean, col_name = 'PBS', fig_name = 'pbs_mean.png')


### Fst Mean
frq = bim.frq
focal = "BRZ"
close = "MES"
outgroup = "EAS"

focal <- frq$frq[CLST == focal, MAF]
close <- frq$frq[CLST == close, MAF]
outgroup <- frq$frq[CLST == outgroup, MAF]

quick_fst <- function(p1, p2) {
  n <- ((p1 - p2)^2) + (((1 - p1) - (1 - p2))^2)
  d <- 2 * ( 1 - ((p1 * p2) + ((1 - p1) * (1 - p2))) )
  res <- n/d
  res[is.nan(res)] <- 0L
  return(res)
}

pbs_fst_mean <- data.table (
  focal.close = sliding_window_mean(quick_fst(focal, close)),
  focal.outgroup = sliding_window_mean(quick_fst(focal, outgroup)),
  close.outgroup = sliding_window_mean(quick_fst(close, outgroup))
)

pbs_fst_mean[, PBS := ( (-log(1L - focal.close)) + (-log(1L - focal.outgroup)) - (-log(1L - close.outgroup)) ) / 2L]
pbs_fst_mean[PBS < 0L, PBS := 0L] # set to 0 negative values

i <- seq(1, (nrow(frq$bim) - 20L), by = 5L) + 19

pbs_fst_mean <- cbind(frq$bim[i, .(CHR, CM, POS, SNP)], pbs_fst_mean[, .(PBS)])
manhattan_plot(pbs_data = pbs_fst_mean, col_name = 'PBS', fig_name = 'pbs_fst_mean.png')
