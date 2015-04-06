# plot global AF histogram for a given population
caf_profiles <- function(pop, breaks, df, hist=FALSE, pop_spec=FALSE) {
  pop_caf <- df[,paste("caf",pop,sep="_")]
  glob_caf <- df[,"caf_all"]
  num_bins <- length(breaks)-1
  pdf <- matrix(,nrow=num_bins)
  bins <- cut(log10(glob_caf), breaks, labels=1:num_bins)
  for (i in 1:(num_bins)){
    bin <- bins == i
    pdf[i] <- sum(pop_caf[bin], na.rm=TRUE)
  }
  return(pdf)
}

# logged <- log10(caf[!colnames(caf)=='gene']+1e-8)
# breaks=seq(from=-7,to=0,by=0.5)
# matplot(prof, type='l',col=pop_colors, lwd=3, lty=1, main="sum POP_CAF over genes in global CAF bin (downsampled)",xaxt='n', xlab='Global CAF', ylab='sum(pop_caf)', bty='n',xaxt='n')
# axis(1, at = 1:length(breaks), labels=breaks, cex.axis = 0.6)
# legend("topleft",pops,fill=pop_colors)


