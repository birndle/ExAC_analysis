AN_CUTOFF = 10000
POP_AN_CUTOFF = 1000
AF_CUTOFF = 0.95
LOW_AF_CUTOFF = 1e-12

# My implementation filters out variant completely if any one population
# fails to satisfy cutoffs
get_caf <- function(df_slice, konrad=TRUE, liberal=F) {
  if (liberal) {
    AN_CUTOFF <- 0
    POP_AN_CUTOFF <- 0
    AF_CUTOFF <- 1
    LOW_AF_CUTOFF <- 0
  }
  
  gene <- df_slice[1,'gene']
  pop_acs <- df_slice[,ac_pop_cols]
  pop_ans <- df_slice[,an_pop_cols]
  ans <- rowSums(pop_ans)
  af <- rowSums(pop_acs)/ans
  pop_af <- pop_acs/pop_ans
  af_df <- cbind(pop_af,af) # af_df has all pops + global af
  colnames(af_df) <- c(paste("caf",pops,sep="_"),"caf_all")
  
  # FILTER VARIANTS AND AGGREGATE ALLELE FREQUENCIES
  cond <- ans > AN_CUTOFF & af > LOW_AF_CUTOFF & af < AF_CUTOFF
  af_df <- af_df[cond,] # Filter variants failing global cutoffs
  pop_ans <- pop_ans[cond,]
  pop_af <- pop_af[cond,]
  if (sum(cond)==0) { # return zero vector if all variants fail global cutoffs
    m <- matrix(0,ncol=ncol(af_df))
    m <- data.frame(m)
    colnames(m) <- c(paste("caf",pops,sep="_"),"caf_all")
    m['gene'] <- gene
    return(m)
  }
  if (!konrad) {  
    # only use variants that pass all population filters
    cond <- ans > AN_CUTOFF & af > LOW_AF_CUTOFF & af < AF_CUTOFF
    cond <- cond & apply(pop_ans, 1, function(row) {all(row > POP_AN_CUTOFF)})
    cond <- cond & apply(pop_af, 1, function(row) {all(row > LOW_AF_CUTOFF & row < AF_CUTOFF)})
    
    af_df <- af_df[cond,]
    caf_df <- data.frame(apply(af_df, 2, sum))
    num_lost <- sum(!cond)
    num_kept <- sum(cond)
  }  else {    
    # Konrad's filter implementation
    # implemented with ac_all and an_all columns
    cond <- pop_ans > POP_AN_CUTOFF & pop_af > LOW_AF_CUTOFF & pop_af < AF_CUTOFF
    if (nrow(af_df)==1){ # only one variant passing
      cond <- c(cond,TRUE) # append a TRUE to the conditions vector for global column
      caf_df <- data.frame(af_df*cond) # zero out pops not passing filters
    } else {
      cond <- cbind(cond,rep(TRUE,nrow(cond))) # append a true column for global
      caf_df <- data.frame(t(colSums(af_df*cond,na.rm=T)))
    }
    caf_df['gene'] <- gene # add gene field 
  }
  return(caf_df)
}

# Split and process multi-gene LoFs 
account_for_multigene <- function(caf_df){
  no_multi_df <- caf_df[grep(",", caf_df$gene, invert=TRUE),] # single gene rows in CAF table
  multi_df <- caf_df[grep(",", caf_df$gene, invert=FALSE),] # multi gene rows
  caf_cols <- colnames(caf_df)[!colnames(caf_df)=='gene']
  for (i in 1:nrow(multi_df)) {
    row <- multi_df[i,]
    # print(class(row$gene))
    genes <- strsplit(as.character(row$gene),",")[[1]]
    for (i in 1:length(genes)) {
      addition <- row[,caf_cols]
      gene_loc_in_df <- no_multi_df$gene==genes[i] # get row position of gene in single-gene CAF table (if present)     
      if (any(gene_loc_in_df)) { # if present in single-gene table
        no_multi_df[gene_loc_in_df, caf_cols] <- no_multi_df[gene_loc_in_df, caf_cols] + addition
      }
      else { # else insert it into table
        addition['gene'] <- genes[i]
        no_multi_df <- rbind(no_multi_df,addition)
      }
    }
  }
  return(no_multi_df)
}

check_against <- function(i){
  gene <- true_caf$Gene[i]
  print(gene)
  k_cols <- paste(toupper(pops),"CAF",sep="_")
  d_cols <- paste("caf",pops,sep="_")
  print(true_caf[true_caf$Gene==gene,k_cols])
  print(d[d$gene==gene,d_cols])
}

# args <- commandArgs(trailingOnly = TRUE)
# if (args[1] == "cluster") {
#   setwd("/humgen/atgu1/fs03/birnbaum/ExAC_analysis")
# } else {
#   setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis")
# }
# source('Rscripts/exac_constants.R')
print("Loading in ExAC ...")
source('Rscripts/load_exac.R')
ac_pop_cols <- paste("ac", pops, sep="_")
an_pop_cols <- paste("an", pops, sep="_")
exac_lof <- subset(exac_all, lof == 'HC' & !(lof_flags %in% c('PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF')))
print("Done.")
library(plyr)

# true_caf <- load_caf_data() # Konrad's CAF
print("Grouping and aggregating ...")
caf <- ddply(exac_lof, 'gene', function(x) { get_caf(x,konrad=TRUE,liberal=T) })
print("Addressing multi-gene variants ...")
caf <- account_for_multigene(caf)
save(caf,file='R_data/caf_lof.Rdata')
# they don't match ... why....
# also deal with variants in multiple genes, i.e. "ENSG10094,ENSG023904"

# lof_dfs <- list.files("downsampled/ac_lof_per_variant",full.names=TRUE)
# load(lof_dfs[1])
# df <- cbind(df, gene=exac_lof$gene) # attach gene column to subsample df
# print('Starting ..')
# caf <- ddply(df, 'gene', function(x) { get_caf(x,konrad=TRUE) })
# caf <- account_for_multigene(caf)
# save(caf, file='Rscripts/caf_test.RData')

# summand_cols <- colnames(caf)[!colnames(caf)=='gene']
# for (i in 2:length(downsampled)) {
#   print(i)
#   load(downsampled[i])
#   df <- cbind(df, gene=exac_lof$gene)
#   new_caf <- ddply(df, 'gene', function(x) {get_caf(x, konrad=TRUE)})
#   new_caf <- account_for_multigene(new_caf)
#   print(all(new_caf$gene == caf$gene))
#   if (!all(new_caf$gene == caf$gene)) {
#     stop("gene order not preserved between samples")
#   }
#   caf[,summand_cols] <- caf[,summand_cols] + new_caf[,summand_cols]
# }
# caf[,summand_cols] <- caf[,summand_cols]/length(downsampled)
# save(caf, file="CAF/downsampled.RData")

# TAKES 13G OF MEMORY TO RUN

# Discrepancies:
# 17310 unique genes in exac_lof vs. 18328 in true_caf
# e.g. ENSG00000183130 not in exac_lof, but is in true_caf





