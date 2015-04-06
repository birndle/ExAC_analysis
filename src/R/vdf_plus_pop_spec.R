vdf <- function(pop, data=exac, breaks=c(-Inf, 0.0005, 0.001, 0.005, 0.01, 0.05, Inf), prob=FALSE, pop_spec=FALSE, prop=FALSE, hom=FALSE, sing=FALSE, doub=FALSE) {
  pop_acs <- data[,paste("ac",pop,sep="_")] # n x i
  pop_ans <- data[,paste("an",pop,sep="_")]
  if (hom){
    pop_homs <- data[,paste("hom",pop,sep="_")]
  }
  pop_af <- pop_acs/pop_ans
  num_bins <- length(breaks)-1
  bins <- cut(pop_af, breaks, labels=1:num_bins)
  if (sing){
    num_bins <- num_bins+1
  } 
  if (doub){
    num_bins <- num_bins+1
  }
  pdf <- matrix(,nrow=num_bins)
  if (pop_spec | hom){
    spec_pdf <- matrix(,nrow=num_bins)
  }
  in_pop <- pop_acs > 0
  unique <- pop_acs == data$ac_adj
  
  # populate singleton and doubleton bins
  if (sing){
    is_singleton <- pop_acs == 1
    pdf[1] <- sum(is_singleton)
    if (doub) { 
      is_doubleton <- pop_acs == 2
      pdf[2] <- sum(is_doubleton)
    }
    if (pop_spec) {
      spec_pdf[1] <- sum(is_singleton & unique)
      if (doub) { spec_pdf[2] <- sum(is_doubleton & unique) }
    }
    else if (hom){
      info1 <- pop_homs[is_singleton]
      if (doub) { info2 <- pop_homs[is_doubleton] }
      if (prop){
        info1 <- info1 > 0
        if (doub) { info2 <- info2 > 0 }
      }
      spec_pdf[1] <- sum(info1)
      if (doub) { spec_pdf[2] <- sum(info2) }
    }
  }
  for (i in 1:(length(breaks)-1)){
    bin <- bins == i & in_pop
    if (sing) {
      i <- i+1
      bin <- bin & !is_singleton
      if (doub) {
        i <- i+1
        bin <- bin & !is_doubleton
      }
    }
    if (pop_spec) {
      spec_pdf[i] <- sum(bin & unique)
    }
    if (hom){
      info <- pop_homs[bin]
      if (prop){
        info <- info > 0
      }
      spec_pdf[i] <- sum(info)
    }
    pdf[i] <- sum(bin)
  }
  if (pop_spec) {
    if (prop) {
      spec_pdf <- spec_pdf / pdf
    }
    if (prob){
      spec_pdf <- spec_pdf/sum(spec_pdf)
    }
    return(spec_pdf)
  }
  else if (hom){
    return(spec_pdf/pdf)
  }
  else {
    if (prob) {
    pdf <- pdf/sum(pdf) 
    }
    return(pdf)
  }
}

plot_vdf <- function(prob=FALSE, norm_by_pop_size=FALSE, log=FALSE, hom=FALSE, prop=FALSE, sing=FALSE) {
  xlab <- c("<.0005", ".0005-.001", ".001-.005", ".005-.01", ".01-.05", ".05+")
  main <- 'Distribution of LoF variants'
  legloc <- "topright"
  if (prob){
    ylab <- "Probability"
    main <- "PDF for LoF variants"
  }
  else if (hom) {
    if (prop) {
      main <- "Proportion of LoF variants with >0 homozygotes"
      ylab <- "Proportion"
      legloc <- "topleft"
    }
    else {
      ylab <- "# homozygous individuals / # variants"
      main <- "Number of homozygotes per LoF variant"
      legloc <- "topleft"
    }
  }
  else {
    ylab <- "Number of LoF variants"
  }
  m <- sapply(pops, function(pop) {vdf(pop, prob=prob, hom=hom, prop=prop, sing=sing)})
  if (norm_by_pop_size) {
    m <- sapply(pops, function(pop) {return(m[,pop]/pop_sizes[pop])})
    main <- "Number of LoFs, normalized by Population size"
  }
  if (log) {
    m <- log10(m)
    ylab <- 'log10(number of LoF variants)'
  }
  matplot(m, type='l', main=main, xlab="af_pop", xaxt="n", ylab=ylab, col=pop_colors, bty='n', lwd=3, lty=1)
  axis(1, at = 1:6, labels=xlab, cex.axis = 0.7)
  legend(legloc, pops, cex=0.8, fill=pop_colors)
}

get_singletons <- function(pop){
  return(sum(exac$ac_adj == 1 & exac[,paste("ac",pop,sep="_")] == 1))
}

plot_pop_specs <- function(log=FALSE, prop=FALSE, prob=FALSE) {
  xlab <- c("<.0005", ".0005-.001", ".001-.005", ".005-.01", ".01-.05", ".05+")  
  m <- sapply(pops, function(pop) {return(vdf(pop, prop=prop, pop_spec=TRUE,prob=prob))})
  ylab <- "# of Variants"
  main <- "Distribution of population-specific LoFs"
  if (log){
    m <- apply(m, 2, log10)
    ylab <- "log10(# of Variants)"
  }
  if (prop){
    main <- "Proportion of LoFs that are population-specific"
    ylab <- "Proportion"
  }
  if (prob){
    main <- "Probability distribution of population-specific LoFs"
    ylab <- "Probability"
  }
  pop_colors <- sapply(pops, function(pop) {return(colors[[pop]])})
  matplot(m, type='l', main=main, xlab="af_pop", xaxt="n", ylab=ylab, col=pop_colors, bty='n', lwd=3, lty=1)
  axis(1, at = 1:6, labels=xlab, cex.axis = 0.7)
  legend("topright",pops, fill=pop_colors, cex=0.8)
}


