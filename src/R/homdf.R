homdf <- function(pop, df=exac) {
  pop_acs <- df[,paste("ac",pop,sep="_")] # n x i
  pop_ans <- df[,paste("an",pop,sep="_")]
  pop_sizes <- pop_ans/2
  pop_homs <- df[,paste("hom",pop,sep="_")]
  hom_freq <- pop_homs/pop_sizes
  breaks <- c(-Inf, 0.0005, 0.001, 0.005, 0.01, 0.05, Inf)
  num_bins <- length(breaks)-1
  pdf <- matrix(,nrow=num_bins)
  bins <- cut(hom_freq, breaks, labels=1:num_bins)
  has_hom <- pop_homs > 0
  for (i in 1:num_bins){
    bin <- bins == i & has_hom
    pdf[i] <- sum(bin)
  }
  return(pdf)
}
setwd('/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis')
