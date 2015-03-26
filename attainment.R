# to run:
# cd /humgen/atgu1/fs03/birnbaum/ExAC_analysis/attainment/
# sh bsub_pops.sh gene <rare/all> 1
# sh bsub_pops.sh lof nope 1
# sh bsub_all.sh uniform gene <rare/all> 1
# sh bsub_all.sh uniform lof nope 1
# sh bsub_all.sh exac gene <rare/all> 1
# sh bsub_all.sh exac lof nope 1

# testing:
# bsub -R "rusage[mem=6]" -q hour -J afr_gene_test -o logs/afr_gene_test Rscript /humgen/atgu1/fs03/birnbaum/ExAC_analysis/Rscripts/attainment.R afr gene /humgen/atgu1/fs03/birnbaum/ExAC_analysis/attainment/curves/afr_gene_test.RData
# bsub -R "rusage[mem=6]" -q hour -J afr_lof_test -o logs/afr_lof_test Rscript /humgen/atgu1/fs03/birnbaum/ExAC_analysis/Rscripts/attainment.R afr lof /humgen/atgu1/fs03/birnbaum/ExAC_analysis/attainment/curves/afr_lof_test.RData 

downsample <- function(pop, max_an, df, prefix) {
  acs <- df[,paste(prefix, pop, sep="_")]
  ans <- df[,paste("an", pop, sep="_")]
  new_ans <- ans*(ans < max_an) + max_an*(ans >= max_an)
  # draw `num_trials` times from a jar of `n` marbles, `true_ac` of which are red 
  take_sample <- function(n, num_trials, true_ac) { 
    return(sum(sample.int(n,num_trials,replace=FALSE) <= true_ac))
  }
  new_acs <- mapply(take_sample, ans, new_ans, acs)
  return(new_acs)
}

an_downsample <- function(pop, max_an, df) {
  ans <- df[,paste("an", pop, sep="_")]
  new_ans <- ans*(ans < max_an) + max_an*(ans >= max_an)
  return(new_ans)
}

# `who` is a vector of pops, i.e. c("afr", "nfe") or just "afr"
lof_attainment <- function(who, sample_sizes, mode='uniform', gene=F, hom=F, rare=F) {
  if (gene) {
    # source("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/Rscripts/gene_lof_burden.R")
    source('/humgen/atgu1/fs03/birnbaum/ExAC_analysis/Rscripts/gene_lof_burden.R')
  }

  check <- function(i) {
    pop <- who[i]
    sample_size <- max_ans[i]
    return(!all(exac_lof[,paste("an",pop,sep="_")] < sample_size))
  }
  profile <- matrix(0,ncol=length(sample_sizes))
  idx <- 1
  for (an in sample_sizes) {
    print(an)
    if (mode == 'uniform') {
      max_ans <- rep.int(round(an/length(who)), length(who))
    } 
    else if (mode == 'exac') {
      sizes <- pop_sizes[who]
      max_ans <- round(sizes/sum(sizes)*an)
    }
    # check that we're not sampling out of the range of any of our pops    
    if (!(all(sapply(1:length(who), check)))) {
      return(profile)
    } 
    print(sprintf("Downsampling to sample size %s ...", an))
    ac_df <- do.call("cbind", mapply(downsample, who, max_ans, MoreArgs = list(df=exac_lof, prefix="ac"), SIMPLIFY=F))
    colnames(ac_df) <- paste("ac", who, sep="_")
    if (hom) {
      hom_df <- do.call("cbind", mapply(downsample, who, max_ans, MoreArgs = list(df=exac_lof, prefix="hom"), SIMPLIFY=F))
      colnames(hom_df) <- paste("hom", who, sep="_")
    }
    an_df <- do.call("cbind", mapply(an_downsample, who, max_ans, MoreArgs = list(df=exac_lof), SIMPLIFY=F))
    colnames(an_df) <- paste("an", who, sep="_")
    print("Done.")

    if (gene) {
      new_df <- data.frame(cbind(ac_df, an_df))
      if (hom) { new_df <- cbind(new_df, hom_df) } # add hom columns
      new_df$gene <- exac_lof$gene # add gene column
      gene_df <- group_by_gene(new_df, hom=hom)
      if (hom) {
        if (rare) {
          col <- "hom_rare"
        } else {
          col <- "hom_all"
        }
      } else if (rare) {
        col <- "ac_rare"
      } else {
        col <- "ac_all"
      }
      profile[idx] <- sum(gene_df[,col] > 0) # number of genes with non-zero entry in given column
    }
    else {
      afs <- rowSums(ac_df)/rowSums(an_df)
      is_rare <- (afs > 0 & afs <= 0.001)
      if (hom) {
        homs <- rowSums(hom_df)
        profile[idx] <- sum(homs > 0 & is_rare, na.rm=T)
      }
      else if (rare) {
        profile[idx] <- sum(is_rare, na.rm=T) # how many rare variants in sample
      }
      else {
        profile[idx] <- sum(afs > 0, na.rm=T)
      }
    }
    print(profile)
  idx <- idx + 1
  }
  return(profile)
}

# e.g. bsub -R "rusage[mem=15]" -q week -J unif -oo attainment/unif_log Rscript Rscripts/attainment.R uniform nope out.RData
args <- commandArgs(trailingOnly = TRUE)
# args[1] == mode of chromsome distribution ("exac" or "uniform" or name of pop)
# args[2] == enter "gene" to indicate that you want to get gene attainment curves
# args[3] == either "rare" or "hom"
  # enter "rare" to indicate you want your attainment curves to only count rare LoFs
    # always enter "rare" for lofs since thats the only analysis we do with them
  # enter "hom" for hom attainment (either gene or variant)
  # enter "all" for gene if you don't want to restrict to "rare" LoFs
# args[4] == output file

if (!is.na(args[1])) {
  setwd("/humgen/atgu1/fs03/birnbaum/ExAC_analysis")
#   setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis")
  print("Loading in ExAC ...")
  source("Rscripts/load_exac.R") # loads exac_all, exac_lof, and pops
  pop_sizes <- apply(exac_lof[,paste("an",pops,sep="_")], 2, max)/2 # in people number
  names(pop_sizes) <- pops
  sample_sizes <- seq(from=1000,to=130000,by=500) # in chromosome number
  by_gene <- args[2] == 'gene'
  only_rare <- args[3] == 'rare'
  hom <- args[3] == 'hom'
  if (args[1] == "uniform") {
    unif_attn_curve <- lof_attainment(pops, sample_sizes, mode='uniform', gene=by_gene, rare=only_rare, hom=hom)
    save(unif_attn_curve, file=args[4])
  } else if (args[1] == "exac") {
    non_unif_attn_curve <- lof_attainment(pops, sample_sizes, mode='exac', gene=by_gene, rare=only_rare, hom=hom)
    save(non_unif_attn_curve, file=args[4])
  } else {
    pop_attn_curve <- lof_attainment(c(args[1]), sample_sizes, gene=by_gene, rare=only_rare, hom=hom)
    save(pop_attn_curve, file=args[4])
  }
}

#########
# OPTIONS
#########

# if (gene)
  # if (hom)
    # if (rare)
      # number of genes with at least one rare, homozyogous LoF - not currently implemented
    # else
      # number of genes with at least one homozygous LoF - testing
  # else
    # if (rare)
      # number of genes with at least one rare LoF - testing
    # else
      # number of genes with at least one LoF - testing
# else 
  # if (hom)
    # number of rare LoFs with at least one homozygous individual - testing
  # if (rare)
    # number of rare LoFs - good
  # else
    # number of LoFs
















