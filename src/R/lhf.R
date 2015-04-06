# load RData containing lhc and hn matrices
# correlate LHF with CAF
# look for genes with significantly high LHF
# stratify by population, then look for genes with significantly high LHF

# setwd('/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/LHC')
# setwd('/humgen/atgu1/fs03/birnbaum/ExAC_analysis/LHC')

get_pop_lhc <- function(lhc) {
  load('R_objects/lhc.RData')
  pops <- c("afr","amr", "eas", "fin", "nfe", "sas")
  pop_table <- read.delim('/humgen/atgu1/fs03/lek/resources/ExAC/ExAC.r0.3_pop_sex.tsv', header=FALSE, row.names=1)
  indivs <- colnames(lhc)
  genes <- rownames(lhc)
  pop_field <- pop_table[indivs,1]
  lhc <- as.matrix(lhc)
  pop_truth <- do.call('cbind', lapply(toupper(pops), function(pop) { pop_field == pop } ))
  colnames(pop_truth) <- pops
  # check that pop_truth has the correct shape and dimensions
  pop_lhc <- lhc %*% pop_truth
  return(pop_lhc)
}

get_lhf <- function(lhc, hn) {
  lhc <- as.matrix(lhc)
  hn <- as.matrix(hn)
  n <- ncol(lhc)
  global_lhc <- lhc %*% rep(1,n)
  global_hn <- hn %*% rep(1,n)
  return(global_lhc/global_hn)
}




