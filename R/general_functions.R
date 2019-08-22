plink <- function(...) {
  values <- c(...)
  args <- names(values)
  plink <- plink_version()
  program_on_path("plink")
  
  if(is.null(args)) args <- ""
  commands <- ifelse(args == "", values, paste(args, values))
  commands <- paste(commands, collapse = " ")
  
  system(paste(plink, commands))
}

plink_version <- function(){
  
  plink_v <- gsub(".*v(\\S..).*", "\\1", system('plink --version', T))
  
  if (plink_v >= 1.9) {
    plink_v <- "plink"
  } else {
    plink_v <- "plink --noweb"
  }
  return(plink_v)
}

read.bim <- function(bim_file, header = F, col.names = c("CHR","SNP","GD","POS","A1", "A2"), colClasses = c("integer", "character", "numeric", "integer", rep("character", 2)), ...) {
  data.table::fread(bim_file, header = F, col.names = col.names, colClasses = colClasses, ...)
}

