# iterate over down-samples
setwd("/humgen/atgu1/fs03/birnbaum/ExAC_analysis")
setwd('/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis')
source('Rscripts/load_exac.R')
ac_dfs <- list.files("downsampled/ac_lof_per_variant",full.names=TRUE)
hom_dfs <- list.files("downsampled/hom_lof",full.names=TRUE)
average_samples <- function(fun, breaks, sing=F, doub=F, prop=F, pdf=F, hom=F, hist=F, pop_spec=F, caf=F, use_median=F, use_mean=F, non_pop_spec=F) {
  num_bins <- length(breaks)-1
  if (sing){
    num_bins <- num_bins + 1
  }
  avg <- matrix(0,nrow=num_bins,ncol=length(pops))
  for (i in 1:length(ac_dfs)) {
    print(i)
    load(ac_dfs[i])
    if (hom) {
      ac_df <- df
      load(hom_dfs[i])
      hom_df <- df
      df <- cbind(ac_df, hom_df[,hom_pop_cols])
    }
    avg <- avg+sapply(pops, function(pop) {return(fun(pop,breaks,df=df,hist=hist,pop_spec=pop_spec,caf=caf,use_median=use_median,use_mean=use_mean,pdf=pdf,non_pop_spec=non_pop_spec))})
  }
  avg <- avg/length(ac_dfs)
  return(avg)
}

pops <- c('amr','nfe','afr','sas','eas','fin')
pop_colors <- sapply(pops, function(pop) {return(colors[[pop]])})
ac_pop_cols <- paste("ac",pops,sep="_")
an_pop_cols <- paste("an",pops,sep="_")
breaks=c(-Inf, 0.001, 0.01, 0.05, Inf)
xlab <- c("Singletons", "<0.001", ".001-.01", ".01-.05", ".05+")
plot_avg <- function(breaks, xlab, pdf=F, prop=F, sing=F, hom=F, log=F) {
  num_bins <- length(xlab) 
  avg <- average_samples(mean_pop_profile, breaks, use_mean=T, sing=T)
  if (log){
    main <- "Distribution of LoFs (downsampled populations)"
    ylab <- "log10(number of variants)"
    avg <- log10(avg)
  }
  ylab <- "mean(POP_AF for LoFs in bin)"
  matplot(log10(avg), type='l', col=pop_colors, main="Mean(POP_AF for LoFs in bin) for LoFs binned by global AF (downsampled)", xlab="global AF", xaxt="n", ylab=ylab, bty='n', lwd=3, lty=1)
  axis(1, at = 1:num_bins, labels=xlab, cex.axis = 0.75)
  legend("topright",pops, fill=pop_colors, cex=0.75)
}




