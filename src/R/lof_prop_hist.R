lof_prop <- function(pop, data) {
  pop_acs <- data[,paste("ac",pop,sep="_")] # n x i
  pop_ans <- data[,paste("an",pop,sep="_")]
  pop_af <- pop_acs/pop_ans
  breaks <- c(-Inf, 0.0005, 0.001, 0.005, 0.01, 0.05, Inf)
  num_bins <- length(breaks)-1
  spectrum <- matrix(,nrow=num_bins)
  bins <- cut(pop_af, breaks, labels=1:num_bins)
  in_pop <- pop_acs > 0
  is_lof <- data$lof == "HC" & is.na(data$lof_flags)
  for (i in 1:num_bins){
    bin <- bins == i & in_pop
    num_vars <- sum(bin)
    num_lofs <- sum(bin & is_lof, na.rm=TRUE)
    spectrum[i] <- num_lofs/num_vars
  }
  return(spectrum)
}
m <- sapply(pops, function(pop) {return(lof_prop(pop, exac_all))})
matplot(m, col=pop_colors,main="Proportion of variants that are LoF",bty='n', lwd=3, xaxt='n',ylab='Proportion', type='l', lty=1, xlab='af_pop')
xlab <- c("<.0005", ".0005-.001", ".001-.005", ".005-.01", ".01-.05", ".05+")
axis(1, at = 1:6, labels=xlab, cex.axis = 0.7)
legend("topright",pops,fill=pop_colors, cex=0.8)


