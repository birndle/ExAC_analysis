# Analyses comparing attainment of LoF gene/variants with respect to size and diversity of sample
library(plyr)

group_by_gene <- function(df, hom=F) {
  grouped <- ddply(df, 'gene', get_burden, hom)
  grouped <- account_for_multigene(grouped, prefix)
  return(grouped)
}

# aggregate LoF variants within a gene slice
get_burden <- function(slice, hom) {
  gene <- slice$gene[1]
  ac_cols <- colnames(slice)[grep("ac", colnames(slice))] # get ac column names
  an_cols <- colnames(slice)[grep("an", colnames(slice))]
  hom_cols <- colnames(slice)[grep("hom", colnames(slice))]
  ac <- rowSums(as.matrix(slice[,ac_cols]))
  an <- rowSums(as.matrix(slice[,an_cols]))
  af <- ac/an
  is_rare <- af <= 0.001 & af > 0 
  df <- data.frame(ac_all=sum(ac), ac_rare=sum(ac*is_rare, na.rm=T), gene=gene)
  if (hom) {
    homs <- rowSums(as.matrix(slice[,hom_cols]))
    df$hom_all=sum(homs)
    df$hom_rare=sum(hom*is_rare, na.rm=T)
  }
  return(df) # df has three fields: "ac_all", "ac_rare", (or "hom"), and "gene"
}

account_for_multigene <- function(df, prefix) {
  no_multi_df <- df[grep(",", df$gene, invert=TRUE),] # single gene rows in CAF table
  multi_df <- df[grep(",", df$gene, invert=FALSE),] # multi gene rows
  counts_cols <- colnames(df)[!colnames(df) == 'gene']
  for (i in 1:nrow(multi_df)) {
    row <- multi_df[i,]
    genes <- strsplit(as.character(row$gene),",")[[1]]
    for (j in 1:length(genes)) {
      counts <- row[,counts_cols]
      loc <- no_multi_df$gene==genes[j] # get row position of gene in single-gene CAF table (if present)     
      if (any(loc)) { # if present in single-gene table
        no_multi_df[loc,counts_cols] <- no_multi_df[loc,counts_cols] + counts
      }
      else { # else insert it into table
        new_row <- data.frame(counts, genes[j])
        names(new_row) <- c(counts_cols, "gene")
        no_multi_df <- rbind(no_multi_df, new_row)
      }
    }
  }
  return(no_multi_df)
}




