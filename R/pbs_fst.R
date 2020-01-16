#' Title
#' 
#' Description
#' 
#' @usage pbs_fst(geno.info, focal, close, outgroup, monomorphic_removal = T, method = "reynolds")
#' 
#' @param plink description
#' @param focal description
#' @param close description
#' @param outgroup description
#' @param filter description
#' @param maf description
#' @param nchrobs description
#' @param monomorphic_removal description
#' @param method description
#' 
#' @return description
#' 
#' @import data.table magrittr
#' @export

pbs_fst <- function(plink, focal, close, outgroup, filter = "FID", maf = NULL, nchrobs = NULL, monomorphic_removal = F, method = "rwc") {
  
  # --- CHECKING ---- #
  
  input <- paste0(plink, "_gdstatIn")
  output <- paste0(plink, "_gdstatOut")
  
  fam_file <- paste0(plink, ".fam")
  fam <- data.table::fread(fam_file)
  
  clst_file <- paste0(output, ".clst")
  clst <- fam[, .(V1, V2)] %>% set_names(c("FID", "ID"))
  clst <- clst[get(filter) %in% c(focal, close, outgroup)]
  
  clst[, state := 'NA']
  clst[get(filter) %in% focal, state := 'focal']
  clst[get(filter) %in% close, state := 'close']
  clst[get(filter) %in% outgroup, state := 'outgroup']
  
  # --- CHECKING ----
  
  # If Plink is installed on system path
  program_on_path("plink")  
  
  # If pre-defined output file's name exists
  if(file.exists(output))
    stop("File ", output, " already exists. Please rename or delete it before running this function.")
  
  # Check method input
  if(!method %in% c("rwc", "wc"))
    stop('Only "rwc" (Reynolds & Cockerham 1983) and "wc" (Weir & Cockerham (1984) methods are available.')
  
  # Check filter input
  if(!filter %in% c("FID", "ID"))
    stop('Only "FID" and "ID" options are available for the filter parameter.')
  
  if(filter == "FID") {
    pops <- c(focal, close, outgroup)
    pops_in_fam <- pops %in% fam[, unique(V1)]
    if(!all(pops_in_fam)) {
      stop("FID [", paste(pops[!pops_in_fam], collapse = ", "), "]", " no present on .fam Plink file!")
    }
  } else {
    pops <- c(focal, close, outgroup)
    pops_in_fam <- pops %in% fam[, V2]
    if(!all(pops_in_fam)) {
      warning("Not all IDs on .fam Plink file!")
    }
  }
  
  # Check maf input
  if(! ( is.null(maf) | (is.numeric(maf) && maf >= 0  && maf <= 1) ) ) {
    stop("MAF must be a numeric value from 0 to 1.")
  } else {
    maf <- ifelse(is.null(maf), "", paste("--maf", maf))
  }
  
  ### COMMON FILTERS ####
  
  # Select target-individuals
  fwrite(clst, clst_file, sep = " ", col.names = F)

  # Extracting target-individuals & Filtering MAF
  if (maf == "") {
    plink(`--bfile` = plink, `--keep` = clst_file, "--keep-allele-order --allow-no-sex", "--make-bed", `--out` = input)
  } else {
    plink(`--bfile` = plink, maf, `--keep` = clst_file, "--keep-allele-order --allow-no-sex", "--make-bed", `--out` = input)
  }
  
  # Remove SNPs with too many missing data
  if(!is.null(nchrobs)) {
    # find SNPs with chromosome observations less than nchrobs
    plink(`--bfile` = plink, `--keep` = clst_file, '--allow-no-sex --freq --within', clst_file, `--out` = input)
    frq.strat <- fread(paste0(input, ".frq.strat"))
    snp_list <- frq.strat[NCHROBS >= nchrobs, unique(SNP)]
    writeLines(snp_list, paste0(input, "_SNPs.txt"))
    # remove those SNPs
    plink(`--bfile` = plink, `--keep` = clst_file, `--extract` = paste0(input, "_SNPs.txt"),
          "--keep-allele-order --allow-no-sex", "--make-bed", `--out` = input)
  }
  
  # Calculate allele frequencies
  if (monomorphic_removal | method == "rwc") {
    
    # Calculate allele frequencies
    plink(`--bfile` = input, "--freq", `--within` = clst_file, "--keep-allele-order --allow-no-sex", "--make-bed", `--out` = output)
    
    # Load files
    bim <- read.bim(paste0(output, ".bim"))
    freq <- data.table::fread(paste0(output, ".frq.strat"))
      
    # Verify and update populations if necessary
    # focal <- ifelse(filter == "FID", focal, "POP1")
    # close <- ifelse(filter == "FID", close, "POP2")
    # outgroup <- ifelse(filter == "FID", outgroup, "POP3")
     
    # Table all alleles frequencies
    freq <- data.table::data.table(POP1 = freq[CLST == "focal", MAF],
                                   POP2 = freq[CLST == "close", MAF],
                                   POP3 = freq[CLST == "outgroup", MAF]
    )
  }
  
  ### WEIR & COCKERHAM 1984 ####
  
  if(method == "wc") {
    
    # Pop1 vs Pop2
    clst[, state := 'NA']
    clst[get(filter) %in% focal, state := 'focal']
    clst[get(filter) %in% close, state := 'close']
    
    data.table::fwrite(clst, clst_file, sep = " ", col.names = F)
    plink(`--bfile` = input, '--fst', `--within` = clst_file, "--keep-allele-order --allow-no-sex", `--out` = output)
    pop1_vs_pop2 <- data.table::fread(paste0(output, ".fst"))
    
    # Pop1 vs Pop3
    clst[, state := 'NA']
    clst[get(filter) %in% focal, state := 'focal']
    clst[get(filter) %in% outgroup, state := 'outgroup']
    
    data.table::fwrite(clst, clst_file, sep = " ", col.names = F)
    plink(`--bfile` = input, '--fst', `--within` = clst_file, "--keep-allele-order --allow-no-sex", `--out` = output)
    pop1_vs_pop3 <- data.table::fread(paste0(output, ".fst"), select = "FST")
    
    # Pop2 vs Pop3
    clst[, state := 'NA']
    clst[get(filter) %in% close, state := 'close']
    clst[get(filter) %in% outgroup, state := 'outgroup']
    
    data.table::fwrite(clst, clst_file, sep = " ", col.names = F)
    plink(`--bfile` = input, '--fst', `--within` = clst_file, "--keep-allele-order --allow-no-sex", `--out` = output)
    pop2_vs_pop3 <- data.table::fread(paste0(output, ".fst"), select = "FST")
    
    # Load fst files
    pbs_fst <- data.table::data.table (
      pop1_vs_pop2[, .(CHR, SNP, POS, FST)],
      pop1_vs_pop3,
      pop2_vs_pop3
    )
    
    # Set column names
    colnames(pbs_fst) <- c("CHR", "SNP", "POS", "POP1.POP2.FST", "POP1.POP3.FST", "POP2.POP3.FST")
    
    # Fix Fst values
    pbs_fst[POP1.POP2.FST < 0L | is.na(POP1.POP2.FST),  POP1.POP2.FST := 0L]
    pbs_fst[POP1.POP3.FST < 0L | is.na(POP1.POP3.FST),  POP1.POP3.FST := 0L]
    pbs_fst[POP2.POP3.FST < 0L | is.na(POP2.POP3.FST),  POP2.POP3.FST := 0L]
    
    # Removing monomorphic alleles in at least two populations
    if(monomorphic_removal) {
      
      m0 <- freq[, lapply(.SD, function(.) . == 0L)][, I := .I]
      m0 <- m0[, M := sum(c(POP1, POP2, POP3)), by = I][, M >= 2L]
      
      m1 <- freq[, lapply(.SD, function(.) . == 1L)]
      m1 <- m1[, I := .I][, M := sum(c(POP1, POP2, POP3)), by = I][, M >= 2L]
      
      mono <- !(m0 | m1)
      cat(" ", sum(!mono), "monomorphic alleles were removed!\n")
      
      pbs_fst <- pbs_fst[mono]

    }
    
    # Calculate PBS
    pbs_fst[, PBS := ( (-log(1L - POP1.POP2.FST)) + (-log(1L - POP1.POP3.FST)) - (-log(1L - POP2.POP3.FST)) ) / 2L]
    pbs_fst[PBS < 0 | is.na(PBS), PBS := 0L][]
  }
  
  #### REYNOLDS & COCKERHAM 1983 ####
  
  if(method == "rwc") {
    
    # Removing monomorphic alleles in at least two populations
    if(monomorphic_removal) {
      
      m0 <- freq[, lapply(.SD, function(.) . == 0L)][, I := .I]
      m0 <- m0[, M := sum(c(POP1, POP2, POP3)), by = I][, M >= 2L]
      
      m1 <- freq[, lapply(.SD, function(.) . == 1L)]
      m1 <- m1[, I := .I][, M := sum(c(POP1, POP2, POP3)), by = I][, M >= 2L]
      
      mono <- !(m0 | m1)
      cat(" ", sum(!mono), "monomorphic alleles were removed!\n")
      
      bim <- bim[mono]
      freq <- freq[mono]
    }
    
    # Reynolds' formula implementation
    fst.reynolds <- function(p1, p2) {
      n <- ((p1 - p2)^2) + (((1 - p1) - (1 - p2))^2)
      d <- 2 * ( 1 - ((p1 * p2) + ((1 - p1) * (1 - p2))) )
      res <- n/d
      res[is.nan(res)] <- 0L
      return(res)
    }
    
    # Calculate Reynolds et al. (1983)'s Fst
    freq[, names(freq) := .(fst.reynolds(POP1, POP2), fst.reynolds(POP1, POP3), fst.reynolds(POP2, POP3))]
    colnames(freq) <- c("POP1.POP2.FST","POP1.POP3.FST", "POP2.POP3.FST")
    
    # Calculate PBS
    freq[, PBS := ( (-log(1L - POP1.POP2.FST)) + (-log(1L - POP1.POP3.FST)) - (-log(1L - POP2.POP3.FST)) ) / 2L]
    freq[PBS < 0 | is.na(PBS), PBS := 0] # set to 0 negative values
    
    # Save object
    pbs_fst <- cbind(bim[, .(CHR, SNP, POS)], freq)
    
  }
  
  files <- c(".bed", ".bim", ".fam", ".log", ".nosex", ".frq.strat", ".fst", ".clst", "_SNPs.txt", "*~")
  unlink(c(paste0(input, files), paste0(output, files)))
  
  pbs_fst[is.infinite(PBS), PBS := 0][]
  return(pbs_fst)
  
}
