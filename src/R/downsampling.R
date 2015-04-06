# max_ans will always be less than or equal to ans, i.e. all(max_ans <= ans) == TRUE
downsample <- function(pop, df, prefix, max_ans) {
  print(pop)
  acs <- df[,paste(prefix, pop, sep="_")]
  ans <- df[,paste("an", pop, sep="_")]
  # new_ans <- ans*(ans < max_an) + max_an*(ans >= max_an) # old
  # draw `num_trials` times from a jar of `n` marbles, `true_ac` of which are red 
  take_sample <- function(n, num_trials, true_ac) { 
    return(sum(sample.int(n,num_trials,replace=FALSE) <= true_ac))
  }
  new_acs <- mapply(take_sample, ans, max_ans, acs)
  return(new_acs)
}

downsample_pops <- function(df, hom=FALSE) {
  if (hom){
    prefix <- "hom"
  }
  else {
    prefix <- "ac"
  }
  max_ans <- apply(df[,an_pop_cols], 1, min)
#   max_an <- min(apply(df[,an_pop_cols],2,max)) # old
  new_df <- data.frame(do.call('cbind', lapply(pops, downsample, df, prefix, max_ans)))
  colnames(new_df) <- paste("ac", pops, sep="_")
  new_df$ans <- max_ans
  return(new_df)
}

# MAIN
# Rscript Rscripts/downsampling.R <output_dir> <num> <lof/syn> <hom?>
args <- commandArgs(trailingOnly = TRUE)
if (!is.na(args[1])) {
  setwd("/humgen/atgu1/fs03/birnbaum/ExAC_analysis")
#   setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis")
  load('Rscripts/exac_pass.RData')
  pops <- c("amr",'nfe','afr','sas','eas','fin')
  an_pop_cols <- paste("an",pops,sep="_")
  ac_pop_cols <- paste("ac",pops,sep="_")
  hom_pop_cols <- paste("hom",pops,sep='_')
  var <- args[1]
  downsample_homs <- args[2] == 'hom'
  output_file <- args[3]
  
  if (var == 'lof') {
    df <- subset(exac, lof == 'HC' & is.na(lof_flags) & use)
  } else if (var == 'syn') {
    df <- subset(exac, consequence == "synonymous_variant" & use)
  }
  downsampled <- downsample_pops(df, hom=downsample_homs)
  save(downsampled, file=output_file)
}
# COMMAND LINE
# for i in {1..100}; do bsub -J syn${i} -q hour -W 2:00 -oo r_logs/sample${i} -R "rusage[mem=9]" Rscript Rscripts/downsampling.R syn nope downsampled/syn_most_recent/sample${i} ; done
