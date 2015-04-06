setwd("/Users/birnbaum88/Desktop/Broad/ExAC_analysis")
load('Rscripts/exac_pass.RData')
exac_all <- exac
exac <- subset(exac_all, lof == 'HC' & is.na(lof_flags) & ac_hom > 0)
hom_pop_cols <- grep("hom_(?!oth)", names(exac), perl=TRUE, value=TRUE)
pops <- unlist(lapply(strsplit(hom_pop_cols, "_"), {function (x) x[[2]]}))

# let k be number of variants (i.e. rows)
get_sample <- function(size, excl_pops) {
  incl <- pops[!pops%in%excl_pops]
#   m <- matrix(,nrow=nrow(exac),ncol=length(incl))
  hom_cols <- paste("hom", incl, sep="_")
  an_cols <- paste("an", incl, sep="_")
  n <- exac[,an_cols]/2
  num_hom <- exac[,hom_cols]
  trial <- n*size # k x num_pops, number of times to sample from population (w/o replacement)
  sum(sapply(1:nrow(exac), function(i){
    ns <- n[i,]
    homs <- num_hom[i,]
    trials <- trial[i,]
#     sapply(1:length(incl), function(j) {
#       m[i,j] <- sum(sample(ns[[j]], trials[[j]], replace=FALSE) < homs[[j]])
#     })
    return(any(sapply(1:length(incl), function(j) { # did any of the pop samples have any homs?
      return(any(sample(ns[[j]], trials[[j]], replace=FALSE) <= homs[[j]])) # did we find any homs in this pop's sample?
    })))
  }))
#   return(m)
} 
samp <- get_sample(1, c("nfe"))

