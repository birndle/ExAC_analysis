setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis")
load('Rscripts/exac_pass.RData')
exac_all <- exac
exac <- subset(exac_all, lof == 'HC' & is.na(lof_flags))
ac_pop_cols <- sort(grep("ac_(?!oth)(?!adj)(?!hemi)(?!het)(?!popmax)(?!hom)", names(exac), perl=TRUE, value=TRUE))
an_pop_cols <- sort(grep("an_(?!oth)(?!adj)(?!hemi)(?!het)(?!popmax)(?!hom)", names(exac), perl=TRUE, value=TRUE))
all_pop_af <- exac[,ac_pop_cols]/exac[,an_pop_cols]
pops <- unlist(lapply(strsplit(ac_pop_cols, "_"), {function (x) x[[2]]}))
pop_sizes <- sapply(pops, function(pop) { 
  return(max(exac[,paste("an",pop,sep="_")])/2)
})

fisher_enrichment <- function(pop,t) {
  print("yo")
  rest <- pops[!pops%in%pop]
  pop_acs <- exac[,paste("ac",pop,sep="_")] # n x i
  pop_ans <- exac[,paste("an",pop,sep="_")]
  rest_acs <- as.matrix(exac[,paste("ac",rest,sep="_")]) %*% rep(1,length(rest))
  rest_ans <- as.matrix(exac[,paste("an",rest,sep="_")]) %*% rep(1,length(rest))
  test_params <- cbind(pop_acs,pop_ans-pop_acs,rest_acs,rest_ans-rest_acs)
  pop_af <- pop_acs/pop_ans
  breaks <- c(-Inf, 0.0005, 0.005, 0.01, 0.05, Inf)
  num_bins <- length(breaks)-1
  p_values <- matrix(,nrow=num_bins+2)
  bins <- cut(pop_af, breaks, labels=3:(length(breaks)+1))
  in_pop <- pop_acs > 0
  unique <- pop_acs == exac$ac_adj
  for (i in 1:(num_bins+2)){
    if (i==1){
      is_singleton <- pop_acs == 1
      # run fisher's exact test on every row of test_params matrix
      p <- apply(test_params[is_singleton,], 1, function(row) {
        return(fisher.test(matrix(row,nrow=2), alternative='greater')$p.value) })
      p_values[i] <- mean(p<t)
    }
    else if (i==2){
      is_doubleton <- pop_acs == 2
      p <- apply(test_params[is_doubleton,], 1, function(row) {
        return(fisher.test(matrix(row,nrow=2), alternative='greater')$p.value) })
      p_values[i] <- mean(p<t)
    }
    else {
      bin <- bins == i & !is_singleton & !is_doubleton & in_pop
      p <- apply(test_params[bin,], 1, function(row) {
        return(fisher.test(matrix(row,nrow=2), alternative='greater')$p.value) })
      p_values[i] <- mean(p<t)
    }
  }
  return(p_values)
}

plot_fisher <- function(t) {
  xlab <- c("Singletons", "Doubletons", "<.05%", ".05-.5%", ".5-1%", "1-5%", "5%+")
  main <- sprintf('Proportion of LoFs with p<%s for Fisher Exact Test',t)
  ylab <- 'Proportion'
  m <- sapply(pops, function(pop) {return(fisher_enrichment(pop,t=t))})
  matplot(m, type='l', main=main, xlab="af_pop", xaxt="n", ylab=ylab,bty='n', lwd=3)
  axis(1, at = 1:7, labels=xlab, cex.axis = 0.7)
  legend("topleft", pops, fill=1:7, cex=0.7)
}


