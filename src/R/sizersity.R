# More analyses comparing attainment of LoF genes with respect to size and diversity of sample
# same ideas, better implementation

lhc_sizersity <- function(levels=c(seq(from=0.1,to=0.9,by=0.1),0.95)) {
  setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/LHC")
  load('R_objects/pop_lhc.RData')
  load('R_objects/pop_hn.RData')
  pops <- c("amr",'nfe','afr','sas','eas','fin')
  max_hns <- apply(pop_hn, 2, max)
  actual_nfe <- max_hns['nfe']
  effec_nfe <- median(max_hns)
  max_hns['nfe'] <- effec_nfe
  distr <- max_hns/sum(max_hns)

  diversity <- matrix(0,ncol=length(levels))
  size <- matrix(0,ncol=length(levels))
  
  do_sample <- function(pop, level, nfe=F) {
    fun <- function(x,y,z) {
      return(sum(sample.int(x, y, replace=F) <= z))
    }
    original <- pop_hn[,pop]
    downed <- round(pop_hn[,"nfe"]*level*distr[pop])
    if (nfe) {
      downed <- round(pop_hn[,"nfe"]*level)
    }
    trials <- mapply(min,original,downed)
    new_counts <- mapply(fun, original, trials, pop_lhc[,pop])
    return(new_counts)
  }
  
  idx = 1
  for (i in levels) {
    print(i)
    print("Downsampling non-NFE pops ..")
    ac <- do.call('cbind', lapply(pops, do_sample, i))
    print("Done.")
    print("Downsampling NFEs ..")
    ac_nfe <- do_sample("nfe", i, nfe=T)
    print("Done.")
    size[idx] <- sum(ac_nfe > 0)
    diversity[idx] <- sum(apply(ac > 0, 1, any))
    idx = idx + 1
    print(size)
    print(diversity)
  }
  return(cbind(t(size), t(diversity)))
}

sizersity <- function(levels) {
  setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis")
  load("R_data/pop_sizes.RData")
  pop_sizes['nfe'] <- median(pop_sizes)
  distr <- pop_sizes/sum(pop_sizes)
  pops <- c("amr",'nfe','afr','sas','eas','fin')
  
  do_sample <- function(pop, level, nfe=F) {
    fun <- function(x,y,z) {
      return(sum(sample.int(x, y, replace=F) <= z))
    }
    original <- exac_lof[,paste("an",pop,sep="_")]
    downed <- round(exac_lof[,"an_nfe"]*level*distr[pop])
    if (nfe) {
      downed <- round(exac_lof[,"an_nfe"]*level)
    }
    trials <- mapply(min,original,downed)
    new_counts <- mapply(fun, original, trials, exac_lof[,paste("ac",pop,sep="_")])
    return(new_counts)
  }

  do_an_sample <- function(pop, level) {
    caps <- round(exac_lof[,"an_nfe"]*level*distr[pop])
    ans <- exac_lof[,paste("an",pop,sep="_")]
    return(mapply(min, ans, caps))
  }

  diversity <- matrix(0,ncol=length(levels))
  size <- matrix(0,ncol=length(levels))

  # weight allele distribution in diverse group by relative sizes of sub-populations
  idx = 1
  t = 0.001
  for (i in levels) {
    print(i)
    # before: for pop in pops, do downsample(acs,ans,round(max_an*an_distr[pop])) where max_an == round(max_nfe*i)
    # now: rather than using one size cap for all variants, determine size cap on per-variant basis
    print("Getting diverse allele counts ..")
    ac <- do.call('cbind', lapply(pops, do_sample, i))
    print("Done.")
    print("Getting diverse chromosome counts ..")
    an <- do.call('cbind', lapply(pops, do_an_sample, i))
    print("Done.")
    pop_ac <- rowSums(ac)
    pop_an <- rowSums(an)
    pop_af <- pop_ac/pop_an

    print("Getting NFE-only allele counts ..")
    ac_nfe <- do_sample("nfe", i, nfe=T)
    print("Done.")
    print("Getting NFE-only chromsome counts ..")
    an_nfe <- round(exac_lof[,"an_nfe"]*i)
    print("Done.")
    af_nfe <- ac_nfe/an_nfe

    size[idx] <- sum(af_nfe > 0 & af_nfe < t, na.rm=T)
    diversity[idx] <- sum(pop_af > 0 & pop_af < t, na.rm=T)
    idx = idx + 1
    print(size)
    print(diversity)
  }
  m <- cbind(t(size),t(diversity))
  return(m)
}

# # run sizersity
# levels <- c(seq(from=0.1,to=0.9,by=0.1),0.95,0.97)
# m <- sizersity(levels)
# matplot(round(levels*max(exac_lof$an_nfe)),m,type='l',xlab='Effective sample size',ylab="# of Rare Variants",main="Accumulation of rare LoFs",xaxt='n',bty='n',lwd=3,lty=1)
# axis(1,at=round(levels*max(exac_lof$an_nfe)))
# legend("topleft",c("NFE only","Mixed"),cex=0.9,fill=1:2)









