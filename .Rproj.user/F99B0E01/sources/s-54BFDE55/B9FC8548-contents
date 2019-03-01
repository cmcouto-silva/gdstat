function (snpIDs, batchsize = 250L, rm.loc = F) {

  if (!is.numeric(length(snpIDs)) || length(snpIDs) == 0)
    stop("Length of SNP vector must be equal or greater than 1.")
  if (!all(grepl("^rs", snpIDs)))
    stop("All SNPs must be codified as Reference SNP ID (starting with 'rs').")

  snp_annot_function <- function(snp_ids) {
    url <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                  "db=snp&id=", paste(snp_ids, collapse = ","), "&report=DocSet",
                  "&tool=LD&email=cmcouto.silva@gmail.com")

    annot <- readLines(url)
    indexes <- grep(pattern = "^rs", annot)
    indexes <- append(x = indexes, values = length(annot) + 1L)
    total <- length(indexes) - 1L
    annot.list <- as.list(replicate(total, NULL))

    for (i in 1L:total) annot.list[[i]] <- annot[indexes[i]:(indexes[i + 1L] - 1L)]

    annot.list <- lapply(annot.list, function(DocSet) {
      index <- grep("GENE=", DocSet, fixed = T)
      gene <- unlist(strsplit(DocSet[index], "GENE=", fixed = T))
      gene <- gene[gene != ""]
    })

    annot.list <- lapply(annot.list, function(x) {
      if (is.null(x)) {
        x <- ""
      } else {
        x
      }
    })
    return(do.call(c, annot.list))
  }

  i <- 0L
  annotation <- character()

  while (i < length(snpIDs)) {
    if ((i + batchsize) > length(snpIDs)) {
      j <- i + abs(length(snpIDs) - i)
    } else {
      j <- i + batchsize
    }

    annotation <- append(annotation, snp_annot_function(snpIDs[(i + 1):j]))
    i <- j
  }

  if(rm.loc) {
    genes <- sapply(annotation, function(annot){
      genes <- unlist(strsplit(annot, ','))
      genes <- genes[!grepl("LOC", genes)]
      paste(genes, collapse = ",")
    })

    annotation <- unname(genes)
    return(annotation)

  } else {
    return(annotation)
  }

}
