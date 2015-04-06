ove <- function(pop, data, use_af=TRUE, median=FALSE, rm_pop_spec=FALSE) {
  pop_acs <- data[,paste("ac",pop,sep="_")] # n x i
  pop_ans <- data[,paste("an",pop,sep="_")] # n x i
  breaks <- c(-Inf, 0.0005, 0.005, 0.01, 0.05, Inf)
  num_bins <- length(breaks)-1
  pdf <- matrix(,nrow=num_bins+2)
  pop_af <- pop_acs/pop_ans
  bins <- cut(pop_af, breaks, labels=3:(length(breaks)+1))
  in_pop <- pop_acs > 0
  unique <- pop_acs == data$ac_adj
  for (i in 1:(num_bins+2)){
    if (i==1){
      is_singleton <- pop_acs == 1
      if (rm_pop_spec) {
        is_singleton <- is_singleton & !unique
      }
      if (!use_af) {
        r <- pop_acs[is_singleton]/data$ac_adj[is_singleton]
      }
      else {
        expected <- (pop_ans[is_singleton]/data$an_adj[is_singleton])*data$ac_adj[is_singleton]
        r <- pop_acs[is_singleton]/expected
      }
    }
    else if (i==2){
      is_doubleton <- pop_acs == 2
      if (rm_pop_spec) {
        is_doubleton <- is_doubleton & !unique
      }
      if (!use_af){
        r <- pop_acs[is_doubleton]/data$ac_adj[is_doubleton]
      }
      else {
        expected <- (pop_ans[is_doubleton]/data$an_adj[is_doubleton])*data$ac_adj[is_doubleton]
        r <- pop_acs[is_doubleton]/expected
      }
    }
    else {
      bin <- bins == i & !is_singleton & !is_doubleton & in_pop
      if (rm_pop_spec) {
        bin <- bin & !unique
      }
      if (!use_af) {
        r <- pop_acs[bin]/data$ac_adj[bin]
      }
      else {
        expected <- (pop_ans[bin]/data$an_adj[bin])*data$ac_adj[bin]
        r <- pop_acs[bin]/expected  
      }     
    }
    if (median){
      pdf[i] <- median(r)
    }
    else {
      pdf[i] <- mean(r)
    }
  }
  return(pdf)
}

plot_ove <- function(median=FALSE, use_af=FALSE, rm_pop_spec=FALSE) {
  xlab <- c("Singletons", "Doubletons","<.05%",".05-.5%", ".5-1%", "1-5%", "5%+")
  if (use_af){
    leg_loc <- "topright"
    if (median) {
      main <- 'Median(observed/expected AC) for LoFs'
      ylab <- "Median(observed/expected AC)"
    }
    else {
      main <- 'Mean(observed/expected AC) for LoFs'
      ylab <- "Mean(observed/expected AC)"
    }
  }
  else {
    leg_loc <- "bottomleft"
    if (median) {
      main <- 'Median(ac_pop/ac_adj) for LoFs'
      ylab <- "Median(ac_pop/ac_adj)"
    }
    else {
      main <- 'Mean(ac_pop/ac_adj) for LoFs'
      ylab <- "Mean(ac_pop/ac_adj)"
    }
  }
  if (rm_pop_spec){
    main <- paste(main," excluding population-specific LoFs", sep=",")
    leg_loc <- "topleft"
  }
  m <- sapply(pops, exac, function(pop) {return(ove(pop, median=median, use_af=use_af, rm_pop_spec=rm_pop_spec))})
  matplot(m, type='l', main=main, xlab="af_pop", xaxt="n", ylab=ylab, bty='n', lwd=3)
  axis(1, at = 1:7, labels=xlab, cex.axis = 0.6)
  legend(leg_loc,pops, fill=1:7, cex=0.65)
}

pop_max <- function(pop){
  pop_acs <- exac[,paste("ac",pop,sep="_")] # n x i
  pop_ans <- exac[,paste("an",pop,sep="_")] # n x i
  pop_af <- pop_acs/pop_ans
  breaks <- c(-Inf, 0.0005, 0.005, 0.01, 0.05, Inf)
  num_bins <- length(breaks)-1
  pdf <- matrix(,nrow=num_bins+2)
  bins <- cut(pop_af, breaks, labels=3:(length(breaks)+1))
  in_pop <- pop_acs > 0
  pop_max <- apply(all_pop_af, 1, function(row) {max(row, na.rm=TRUE)})
  for (i in 1:(num_bins+2)){
    if (i==1) {
      is_singleton <- pop_acs == 1
      pdf[i] <- mean(pop_max[is_singleton] == pop_af[is_singleton])
    }
    else if (i==2){
      is_doubleton <- pop_acs == 2
      pdf[i] <- mean(pop_max[is_doubleton] == pop_af[is_doubleton])
    }
    else {
      bin <- bins == i & !is_singleton & !is_doubleton & in_pop
      pdf[i] <- mean(pop_max[bin] == pop_af[bin])
    }
  }
  return(pdf)
}

plot_pop_max <- function() {
  xlab <- c("Singletons", "Doubletons","<.05%",".05-.5%", ".5-1%", "1-5%", "5%+")
  m <- sapply(pops, pop_max)
  pop_colors <- sapply(pops, function(pop) {return(colors[[pop]])})
  matplot(m, type='l', main="Proportion of LoFs for which af_pop == pop_max", xlab="af_pop", xaxt="n", ylab="Proportion", col=pop_colors, bty='n', lwd=3)
  axis(1, at = 1:7, labels=xlab, cex.axis = 0.7)
  legend("bottomleft",pops, fill=pop_colors, cex=0.8)
}
