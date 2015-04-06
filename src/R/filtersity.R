filtersity <- function(use_median=F, use_max=F) {
  library(gtools)
  lof_dfs <- list.files("downsampled/ac_lof_per_variant",full.names=TRUE)
  load(lof_dfs[1])
  mean_df <- df
  print("Averaging over downsampled dfs ..")
  for (i in 2:length(lof_dfs)) {
    load(lof_dfs[i])
    mean_df <- mean_df + df
  }
  mean_df <- mean_df/length(lof_dfs)
  print('Done.')
  # proportion of rare LoFs filtered by pop_max
  get_proportion <- function(pops,df,ub=0.01,lb=0.05) {
    acs <- df[,paste('ac',pops,sep='_')]
    ans <- df[,paste('an',pops,sep='_')]
    afs <- acs/ans
    glob_af <- rowSums(acs)/rowSums(ans)
    cands <- glob_af > 0 & glob_af < ub
    before <- sum(cands,na.rm=T)
    filter <- apply(afs, 1, function(x) { any(x > lb) })
    after <- sum(cands & filter, na.rm=T) # number of candidates that fail filter
    return(after/before) # return proportion of candidates filtered
  }
  num_pops <- length(pops)
  levels <- seq(from=0.1,to=1,by=0.1)
  m <- matrix(0,ncol=length(pops)-1,nrow=length(levels))
  for (i in levels) {
    print(i)
    new_df <- data.frame(matrix(,nrow(mean_df),length(cols),dimnames=list(c(),cols)))
    # downsample data frame to current level
    for (pop in pops) {
      acs <- mean_df[,paste("ac",pop,sep='_')]
      ans <- mean_df[,paste("an",pop,sep='_')]
      l <- downsample(acs,ans,round(max(ans)*i)) # use mean_df for downsampling parameters
      new_df[,paste("ac",pop,sep="_")] <- l[[1]]
      new_df[,paste("an",pop,sep="_")] <- l[[2]]
    }
    for (j in 2:length(pops)) {
      print(j)
      groups <- combinations(num_pops, j, pops)
      p <- apply(groups, 1, function(g) { get_proportion(g,new_df) })
      if (use_median) {
        m[i*10,j-1] <- median(p)
      }
      else if (use_max) {
        m[i*10,j-1] <- max(p)
      } 
      else {
        m[i*10,j-1] <- mean(p)
      }
    }
    print(m)
  }
  return(m)
}

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
use_median <- mode == 'median'
use_max <- mode == 'max'

m <- filtersity(use_median=use_median, use_max=use_max)
save(m, file=args[2])

# # run filtersity
# m <- filtersity(use_max=T)
# levels <- seq(from=0.1,to=1,by=0.1)
# matplot(levels,m,type='l',lty=1,lwd=3,main='Proportion of rare LoFs filtered by pop_max',xlab='Size of sample', ylab='Proportion filtered',xaxt='n',bty='n')
# axis(1,at=levels)
# legend("bottomright",c('2 pop','3 pops', '4 pops', '5 pops', '6 pops'),fill=1:5, cex=0.7, horiz=T) 
