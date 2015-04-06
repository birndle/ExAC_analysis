pop_af_hist <- function(df, pops, down=F, sing=T, breaks=c(-Inf,0.001,0.01,0.05,Inf), pop_spec=F, ub=0) {
  pop_acs <- data.frame(df[,paste("ac",pops,sep="_")])
  if (down) {
    pop_afs <- pop_acs/df$ans
  } else {
    pop_ans <- data.frame(df[,paste("an",pops,sep="_")])
    pop_afs <- data.frame(pop_acs/pop_ans)
  }
  colnames(pop_afs) <- paste("af",pops,sep="_")
  num_bins <- length(breaks)-1 + 1 # +1 for singleton bin
  spectrum <- matrix(,nrow=num_bins,ncol=length(pops))
  bins <- apply(pop_afs, 2, function(afs) {cut(afs, breaks, labels=2:num_bins)})
  unique <- apply(pop_acs, 2, function(x) { x == rowSums(pop_acs) })
  if (ub & pop_spec) {
    return (colSums(unique & (pop_afs < ub), na.rm=TRUE))
  }
  if (sing) {
    start <- 1
  } else {
    start <- 2
  }
  
  for (i in start:num_bins) {
    if (i==1) {
      is_singleton <- pop_acs == 1
      if (pop_spec) {
        is_singleton <- is_singleton & unique
      }
      spectrum[i,] <- colSums(is_singleton, na.rm=T)
    }
    else {
      bin <- bins == i & pop_acs > 0
      if (sing) {
        bin <- bin & !is_singleton
      }
      if (pop_spec) {
        bin <- bin & unique
      }
      spectrum[i,] <- colSums(bin, na.rm=T)
    }
  }
  if (!sing) {
    spectrum <- spectrum[2:nrow(spectrum),]
  }
  return(spectrum)
}

plot_pop_af_hist <- function(df, title, down=F, log=F) {
  source("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/Rscripts/exac_colors.R")
  xticks <- c("Singletons", "< 0.001", "0.001-0.01", "0.01-0.05", "0.05+")
  m <- pop_af_hist(df, pops, down=down)
  ylab <- "Number of variants"
  if (log) { 
    m <- log10(m) 
    ylab <- sprintf("log10(%s)", ylab)
  }
  matplot(m, col=pop_colors, type='l', lty=1, lwd=3, main=title, xlab="population AF", ylab=ylab, xaxt='n', bty='n')
  axis(1, at=1:5,xticks)
  legend("topright", pops, fill=pop_colors)
}

hists1 <- pop_af_hist(exac_all, c("nfe", "fin"), sing=F)
hists2 <- pop_af_hist(exac_lof_all, c("nfe", "fin"), sing=F)
hists3 <- pop_af_hist(exac_lof_use, c("nfe", "fin"), sing=F)
hists4 <- pop_af_hist(exac_lof_use_down, c("nfe", "fin"), sing=F, down=T)
m <- cbind(hists1[,2]/hists1[,1], hists2[,2]/hists2[,1], hists3[,2]/hists3[,1], hists4[,2]/hists4[,1], c(1,1,1,1))
matplot(m, type='l', col=c('red', 'green', 'blue', 'orange', 'black'), lty=c(1,1,1,1,2), lwd=3, xaxt='n', xlab='af_pop', ylab='# FIN variants / # NFE variants', main='FINs vs. NFEs')
xticks <- c("< 0.001", "0.001-0.01", "0.01-0.05", "0.05+")
axis(1, at=1:4, xticks)
legend("bottomright", fill=c('red', 'green', 'blue', 'orange'), c("all variants", "HC LoFs", "HC LoFs with use", "HC LoFs with use, downsampled"), cex=0.7)

lof_filter <- function(df,pop,lb=0.01,ub=0.001,hgmd=F,nfe=F) {
  ac_col <- paste("ac",pop,sep="_")
  an_col <- paste("an",pop,sep="_")
  pop_ac <- df[,ac_col]
  pop_an <- df[,an_col]
  rest_ac <- df[,ac_pop_cols[!ac_pop_cols==ac_col]]
  rest_an <- df[,an_pop_cols[!an_pop_cols==an_col]]
  rest_af <- rowSums(rest_ac)/rowSums(rest_an)
  cands <- rest_af < ub & rest_af > 0
  if (hgmd) {
    cands <- cands & in_hgmd
  }
  before <- sum(cands, na.rm=T)
  filter <- (pop_ac/pop_an) > lb
  after <- sum(cands & filter, na.rm=T) # number of LoFs failing pop_af filter
  if (nfe) {
    nfe_uniq <- rowSums(rest_ac) == df[,'ac_nfe']
    after1 <- sum(cands & filter & nfe_uniq, na.rm=T)
    return(after1/before)
  }
  return(after/before)
}

hom_filter <- function(df,hom_df,pop,ub=0.001,hgmd=F) {
  homs <- hom_df[,paste('hom',pop,sep='_')]
  pop_acs <- df[,ac_pop_cols]
  pop_ans <- df[,an_pop_cols]
  glob_af <- rowSums(pop_acs)/rowSums(pop_ans)
  cands <- rest_af < ub & rest_af > 0
  if (hgmd) {
    cands <- cands & in_hgmd
  }
  before <- sum(cands, na.rm=T)
  filter <- homs > 0
  after <- sum(cands & filter, na.rm=T) # number of LoFs failing pop_af filter
  return(after/before)
}

# load in and average together downsampled df's
load_downsampled_dfs <- function(lof=T, cluster=F) {
  if (cluster) {
    setwd('/humgen/atgu1/fs03/birnbaum/ExAC_analysis/')
  } else {
    setwd('/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/')
  }
  if (lof) {
    dfs <- list.files('downsampled/lof_most_recent', full.names=T)
  } else {
    dfs <- list.files('downsampled/syn_most_recent', full.names=T)
  }
  tot <- get(load(dfs[1]))
  for (i in 2:length(dfs)) {
    add <- get(load(dfs[i]))
    tot <- tot + add  
  }
  tot <- tot/length(dfs)
  return(round(tot))
}

lof_vs_syn <- function(plot=F) {
#   load("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/R_data/exac_pass.RData")
#   exac_syn <- subset(exac, consequence == "synonymous_variant")
#   exac_lof <- subset(exac, lof == 'HC' & is.na(lof_flags))
  source("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/Rscripts/exac_colors.R")
  breaks <- c(-Inf, 0.001, 0.01, 0.05, Inf)
  lof_df <- load_downsampled_dfs(lof=T)
  lof <- pop_af_hist(lof_df, pops, down=T, breaks=breaks)
  syn_df <- load_downsampled_dfs(lof=F)
  syn <- pop_af_hist(syn_df, pops, down=T, breaks=breaks)
  if (plot) {
    xlab <- c("Singletons", "<0.001", ".001-.01", ".01-.05", ".05+")      
    matplot(lof/syn, type='l',col=pop_colors,main="LoF/Synonymous (downsampled)", xlab='pop_af', ylab='# LoF / # Synonymous', xaxt='n', bty='n',lwd=3,lty=1)
    axis(1, at = 1:length(xlab), labels=xlab, cex.axis = 0.75)
    legend("topright",pops, fill=pop_colors, cex=0.75)
  }
  return(lof/syn)
}





