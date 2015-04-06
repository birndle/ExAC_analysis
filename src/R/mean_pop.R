# multi-purpose function for generating histograms with LoFs binned by global AF

mean_pop_profile <- function(pop, breaks, df=exac, glob=T,pdf=F,use_mean=F,use_median=F,hist=F,pop_spec=F,caf=F,non_pop_spec=F) {
  pop_acs <- df[,paste("ac",pop,sep="_")]
  pop_ans <- df[,paste("an",pop,sep="_")] 
  pop_af <- pop_acs/pop_ans
#   glob_ac <- rowSums(df[,ac_pop_cols])
#   glob_an <- rowSums(df[,an_pop_cols])
#   glob_af <- glob_ac/glob_an
  num_bins <- length(breaks)-1 + 1 # +1 for singleton bin
  spectrum <- matrix(,nrow=num_bins)
  bins <- cut(glob_af, breaks, labels=2:num_bins)
  in_pop <- pop_acs > 0
  unique <- pop_acs == glob_ac
  if (use_mean) {
    fun <- mean
  }
  else if (use_median) {
    fun <- median
  }
  else if (caf) {
    fun <- sum
  }
  for (i in 1:(num_bins)){
    if (i==1){
      is_singleton <- glob_ac == 1
      if (hist | pdf){
        spectrum[i] <- sum(is_singleton & in_pop, na.rm=TRUE)
      }
      else if (pop_spec) {
        spectrum[i] <- sum(is_singleton & unique) # for singletons, being unique and being in the pop are synonymous
      }
      else {
        spectrum[i] <- fun(pop_af[is_singleton], na.rm=TRUE) # NAs occur when pop_an == 0        
      }      
    }
    else {
      bin <- bins == i & !is_singleton
      if (hist | pdf) {
        spectrum[i] <- sum(bin & in_pop, na.rm=TRUE)
      }
      else if (pop_spec) {
        spectrum[i] <- sum(bin & unique)
      }
      else if (non_pop_spec) {
        if (caf | use_median | use_mean) {
          spectrum[i] <- fun(pop_af[bin & in_pop & !unique])
        }
        else{
          spectrum[i] <- sum(bin & in_pop & !unique)
        }
      }
      else {
        spectrum[i] <- fun(pop_af[bin], na.rm=TRUE)
      }
    }
  }
  if (pdf) {
    spectrum <- spectrum/sum(spectrum)
  }
  return(spectrum)
}




