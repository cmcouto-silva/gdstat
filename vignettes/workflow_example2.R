#### Main Analysis ####
####################################################################################################!

# Load .bim / .frq.strat files
eur_afr <- fread("~/cmcouto.silva@usp.br/lab_files/all_datasets/1000G/vcf/LCT/eur_vs_afr.fst")
eur_eas <- fread("~/cmcouto.silva@usp.br/lab_files/all_datasets/1000G/vcf/LCT/eur_vs_eas.fst")
afr_eas <- fread("~/cmcouto.silva@usp.br/lab_files/all_datasets/1000G/vcf/LCT/afr_vs_eas.fst")

eur_afr[is.nan(FST) | FST < 0, FST := 0]
eur_eas[is.nan(FST) | FST < 0, FST := 0]
afr_eas[is.nan(FST) | FST < 0, FST := 0]

pbs2 <- ((-log(1 - eur_afr$FST)) + (-log(1 - eur_eas$FST)) - (-log(1 - afr_eas$FST))) / 2
pbs2[pbs2 < 0] <- 0L

pbs2 <- data.table(POS = eur_afr[, POS], EUR_AFR = eur_afr[, FST], EUR_EAS = eur_eas[, FST], 
                  AFR_EAS = afr_eas[, FST], PBS_EUR = pbs2)

pbs2[POS == LCT]

PBS_distrQT2 <- quantile(pbs2[, PBS_EUR], c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
PBS_distrQT2 # estimar os quantils
  
# PBS for individuals SNPs
pbs <- pbs_fst(geno.info = bim.frq, focal = "EUR", close = "EAS", outgroup = "AFR", monomorphic_removal = F)
# manhattan_plot(pbs_data = pbs, col_name = 'PBS', fig_name = 'pbs_whole.png')


reynolds_fst(bim.frq, "AFR", "EAS")[lct_pos]
reynolds_fst(bim.frq, "AFR", "EUR")[lct_pos]
reynolds_fst(bim.frq, "EAS", "EUR")[lct_pos]


# PBS for mean of SNPs
pbs_mean <- sliding_window_mean_dt(pbs, 'PBS')
pbs_mean_peaks <- get_peaks(pbs_mean, quantile = 0.999, snp.annot = F, rm.loc = F)



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