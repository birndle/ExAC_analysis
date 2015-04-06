# Analyses using LoF haplotype counts in ExAC. 

# load in and average together downsampled df's
load_downsampled_dfs <- function(cluster=T) {
  if (cluster) {
    setwd('/humgen/atgu1/fs03/birnbaum/ExAC_analysis/LHC')
  } else {
    setwd('/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/LHC')
  }
  dfs <- list.files('downsampled/second', full.names=T)
  load(dfs[1])
  lhc_df <- m
  for (i in 2:length(dfs)) {
    load(dfs[i])
    lhc_df <- lhc_df + m   
  }
  lhc_df <- lhc_df/length(dfs)
  return(lhc_df)
}

run_compare_combos <- function(plot=F) {
  best <- plot_combos(7,get_best=T)
  print(best)
  lhc_df <- load_downsampled_dfs(cluster=F)
  m <- compare_combos(lhc_df, best)
  if (plot) {
    levels <- seq(from=0.1,to=1,by=0.1)
    matplot(levels,m,type='l',lty=1,lwd=3,xlab='sample size',ylab='# of Genes discovered with LoF',main="Discovery of Genes with LoFs")
    legend("bottomright",best,fill=1:6,cex=0.8)
  }
  return(m)
}

# best is a vector where the i'th element is the best combination of i populations 
# accepts downsampled df's
compare_combos <- function(df, best, levels=seq(from=0.1,to=1,by=0.1)) {  
  load('R_objects/pop_hn.RData')
  size_cap <- min(apply(pop_hn, 2, max))
  profile <- matrix(0,nrow=length(levels),ncol=length(best))
  pops <- c('afr', 'fin', 'nfe', 'amr', 'sas', 'eas')
  
  do_sample <- function(pop, level, k, nfe=F) {
    fun <- function(x,y,z) {
      return(sum(sample.int(x, y, replace=F) <= z))
    }
    original <- pop_hn[,pop]
    original <- original*(original <= size_cap) + size_cap*(original > size_cap) # recall the downsampled HN
    trials <- round(original*level/k)
    new_counts <- mapply(fun, original, trials, df[,toupper(pop)])
    return(new_counts)
  }

  for (i in 1:length(levels)) {
    print(profile)
    for (j in 1:length(best)) {
      print(j)
      use_pops <- strsplit(best[j],",")[[1]]
      print(use_pops)
      ac <- do.call('cbind', lapply(tolower(use_pops), do_sample, levels[i], j))
      profile[i,j] <- sum(apply(ac > 0, 1, any))
    }
  }
  colnames(profile) <- best
  return(profile)
}


# compare gene-LoF discovery for different, k-sized combinations of populations 
# accepts downsampled df's
best_combo <- function(df, k, levels=seq(from=0.1,to=1,by=0.1)) {
  library(gtools)
  load('R_objects/pop_hn.RData')
  size_cap <- min(apply(pop_hn, 2, max))
  pops <- toupper(c('afr', 'fin', 'nfe', 'amr', 'sas', 'eas'))
  combos <- combinations(length(pops), k, pops) 
  
  do_sample <- function(pop, level) {
    fun <- function(x,y,z) {
      return(sum(sample.int(x, y, replace=F) <= z))
    }
    original <- pop_hn[,tolower(pop)]
    original <- original*(original <= size_cap) + size_cap*(original > size_cap) # recall the downsampled HN
    trials <- round(original*level)
    new_counts <- mapply(fun, original, trials, df[,pop])
    return(new_counts)
  }
  
  profile <- matrix(0,nrow=length(levels),ncol=nrow(combos))
  colnames(profile) <- apply(combos,1,paste,collapse=",")
  for (i in 1:length(levels)) {
    print(profile)
    for (j in 1:nrow(combos)) {
      print(j)
      use_pops <- combos[j,]
      ac <- do.call('cbind', lapply(use_pops, do_sample, levels[i]))  
      profile[i,j] <- sum(apply(ac > 0, 1, any))
    }
  }
  return(profile)
}

# plot discovery profiles for all nCr combinations of pops, given r
plot_combos <- function(r, get_best=F) {
  setwd('/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/LHC')
  xax <- seq(from=0.1,to=1,by=0.1)
  xlab <- "sample size"
  ylab <- "number of genes found with LoF"
  main <- "Discovery of LoF genes"
  
  if (get_best) {
    evaluate <- function(t) {
     t <- apply(t, 2, sum) 
     return(names(t)[t==max(t)])
    }
    best <- matrix(0,ncol=6)
    nums <- c("one","two","three","four","five","six")
    for (i in 1:length(nums)) {
      load(sprintf('R_objects/combo_profiles/second/%s.RData', nums[i])) # load in combo
      best[i] <- evaluate(combo)      
    }
    return(best)
  } 
  # where r == 1
  source('../Rscripts/exac_colors.R')
  if (r == 1) {
    load('R_objects/combo_profiles/one.RData')
    clrs <- pop_colors[colnames(x)]
    matplot(xax,x, main=main,type='l',col=clrs,lwd=3,lty=1,xlab=xlab,ylab=ylab)
    legend("topleft",colnames(x),fill=clrs,cex=0.85,ncol=2)
  }
  if (r == 2) {
    # where r == 2
    load('R_objects/combo_profiles/two.RData')
    matplot(xax,x,type='p',col=1:6,pch=1:ncol(x),lty=2,xlab=xlab,ylab=ylab,main=main)
    legend("bottomright",colnames(x),col=1:6,pch=1:ncol(x), cex=0.75, ncol=2)
  }
  if (r == 3) {
    load('R_objects/combo_profiles/three.RData')
    matplot(xax,x,type='p',col=pop_colors,pch=1:ncol(x),lty=2,xlab=xlab,ylab=ylab,main=main)
    legend("bottomright",colnames(x),col=pop_colors,pch=1:ncol(x), cex=0.7, ncol=3)
  }
  if (r == 4) {
    load('R_objects/combo_profiles/four.RData')
    matplot(xax,x, type='p',col=pop_colors,pch=1:ncol(x),lty=2,xlab=xlab,ylab=ylab,main=main)
    legend("bottomright",colnames(x),col=pop_colors,pch=1:ncol(x), cex=0.7, ncol=2)
  }
  if (r == 5) {
    # where r == 5
    load('R_objects/combo_profiles/five.RData')
    matplot(xax,x,type='l',col=pop_colors,lty=1,lwd=3,xlab=xlab,ylab=ylab,main=main)
    legend("bottomright",colnames(x),fill=pop_colors, cex=0.85, ncol=1)
  }
  if (r == 6) {
    # where r == 6
    load('R_objects/combo_profiles/six.RData')
    plot(x,type='l',lwd=3,col='red',xlab=xlab,ylab=ylab,main=main)
    legend("bottomright",'all populations',fill='red')
  }
}
# run main script to generate downsampled dfs
args <- commandArgs(trailingOnly = TRUE)
if (!is.na(args[1])) {
  # downsample pop_lhc tables to FIN size
  # Rscript lhf_analyses.R downsample <i>
  if (args[1] == 'downsample') {
    # run from within LHC directory
    setwd('/humgen/atgu1/fs03/birnbaum/ExAC_analysis/LHC')
    pops <- c('afr', 'fin', 'nfe', 'amr', 'sas', 'eas')

    print("Loading in HN table ...")
    load('R_objects/pop_hn.RData')
    print("Done.")
    print("Loading in LHC table ...")
    load('R_objects/pop_lhc_floored.RData')
    print("Done.")
    
    # downsample lhc tables
    df <- pop_lhc_floored
    level <- min(apply(pop_hn, 2, max))
    # # for-loop implementation takes ~6 min
    # double lapply takes ~5 min
    cols <- colnames(df)[colnames(df)!='OTH']
    do_sample <- function(row, col) {
      pop <- col
      hn <- pop_hn[row,tolower(pop)]
      return(sum(sample.int(hn,min(c(level,hn)),replace=FALSE) <= df[row,col]))
    }
    print("downsampling..")
    m <- do.call('rbind', lapply(rownames(df), function(row) { 
      do.call('cbind', lapply(cols, function(col) { do_sample(row,col) }))}))
    print("Done.")
    rownames(m) <- rownames(df)
    colnames(m) <- cols
    save(m, file=sprintf("downsampled/second/sample%s.RData",args[2]))    
  }
  # compute profile of gene-LoF accumulation for nCk different population subsets  
  # Rscript lhf_analyses.R combos <k> <RData_output_file>
  else if (args[1] == 'combos') {
    print('Loading in downsampled pop LHC matrix ...')
    df <- load_downsampled_dfs(cluster=T)
    print("Done.")
    print(head(df))
    print("Computing profiles ...")
    combo <- best_combo(df, as.double(args[2]), levels=seq(from=0.1,to=1,by=0.1))
    print("Done.")
    # print(combo)
    save(combo, file=args[3])
  }
}

# plot combo_profile
matplot(round(levels*size_cap), combo_profile,type='l',xlab='Sample size in people',ylab="# of Genes with at least one LoF",main="Attainment of genes with LoFs",xaxt='n',bty='n',lwd=3,lty=1)
axis(1,at=round(levels*size_cap))
legend("bottomright",colnames(combo_profile),fill=1:6, cex=0.7)


