setwd("/humgen/atgu1/fs03/birnbaum/ExAC_analysis")
#   setwd("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis")
print("Loading in ExAC ...")
source("Rscripts/load_exac.R") # loads exac_all, exac_lof, and pops
pop_sizes <- apply(exac_lof[,paste("an",pops,sep="_")], 2, max)/2 # in people number
names(pop_sizes) <- pops
sample_sizes <- seq(from=1000,to=121000,by=500) # in chromosome number

# min_ans <- apply(exac_lof[,paste("an", pops, sep="_")], 1, min)

# to run:
# cd /humgen/atgu1/fs03/birnbaum/ExAC_analysis/attainment/
# sh bsub_pops.sh gene <rare/all/hom> 1a
# sh bsub_pops.sh lof <rare/all/hom> 1a
# sh bsub_all.sh uniform gene <rare/all/hom> 1a
# sh bsub_all.sh uniform lof <rare/all/hom> 1a
# sh bsub_all.sh exac gene <rare/all/hom> 1a
# sh bsub_all.sh exac lof <rare/all/hom> 1a

# check directories:
# for d in ./*; do echo $d; for d2 in ${d}/*; do ls $d2; done; done # in curves/pops
# for d in ./*; do ls $d; done # in curves/all/uniform

# new strategy:
# for sample size s and variant v, sample all populations down to min(s, min(pop_ans(v))) 

downsample <- function(pop, max_an, df, prefix, mixed=F) {
  which_pops <- function(level) {
    return(pops[pop_sizes*2 >= level])
  }
  print(mixed)
  acs <- df[,paste(prefix, pop, sep="_")]
  ans <- df[,paste("an", pop, sep="_")]

  if (!mixed) {
    use_pops <- which_pops(max_an)
    print(use_pops)
    min_ans <- apply(as.matrix(exac_lof[,paste("an", use_pops, sep="_")]), 1, min)
    new_ans <- min_ans*(min_ans < max_an) + max_an*(min_ans >= max_an) # new
  }
  else {    
    new_ans <- ans*(ans < max_an) + max_an*(ans >= max_an) # old
  }

  # draw `num_trials` times from a jar of `n` marbles, `true_ac` of which are red 
  take_sample <- function(n, num_trials, true_ac) { 
    return(sum(sample.int(n,num_trials,replace=FALSE) <= true_ac))
  }
  new_acs <- mapply(take_sample, ans, new_ans, acs)
  return(new_acs)
}

an_downsample <- function(pop, max_ans, df) {
  ans <- df[,paste("an", pop, sep="_")]
  new_ans <- ans*(ans < max_ans) + max_ans*(ans >= max_ans)
  return(new_ans)
}

# `who` is a vector of pops, i.e. c("afr", "nfe") or just "afr"
lof_attainment <- function(who, sample_sizes, mode='uniform', gene=F, hom=F, rare=F, get_df=F) {
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
    } else if (mode == 'exac') {
      sizes <- pop_sizes[who]
      max_ans <- round(sizes/sum(sizes)*an)
    } else if (mode == 'world') {
      sizes <- world_pop_sizes[who]
      max_ans <- round(sizes/sum(sizes)*an)
    }  
    mixed <- length(who) != 1

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
        col <- "lof_rare"
      } else {
        col <- "lof_all"
      }
      profile[idx] <- sum(gene_df[,col] > 0) # number of genes with non-zero entry in given column
    }
    else {
      afs <- rowSums(ac_df)/rowSums(an_df) # compute af across all pops in sample
      is_rare <- (afs > 0 & afs <= 0.01)
      if (hom) {
        homs <- rowSums(hom_df)
        profile[idx] <- sum(homs > 0, na.rm=T) # how many homozygous LoF sites
      }
      else if (rare) {
        profile[idx] <- sum(is_rare, na.rm=T) # how many rare LoF sites
      }
      else {
        profile[idx] <- sum(afs > 0, na.rm=T) # how many LoF sites
      }
    }
    print(profile[1:idx])
  idx <- idx + 1
  }
  return(profile)
}

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
  by_gene <- args[2] == 'gene'
  only_rare <- args[3] == 'rare'
  hom <- args[3] == 'hom'
  if (args[3] == 'hom_rare') {
    only_rare <- T
    hom <- T
  }

  if (args[1] == "uniform") {
    unif_attn_curve <- lof_attainment(pops, sample_sizes, mode='uniform', gene=by_gene, rare=only_rare, hom=hom)
    save(unif_attn_curve, file=args[4])
  } else if (args[1] == "exac") {
    non_unif_attn_curve <- lof_attainment(pops, sample_sizes, mode='exac', gene=by_gene, rare=only_rare, hom=hom)
    save(non_unif_attn_curve, file=args[4])
  } else if (args[1] == 'world') {
    world_attn_curve <- lof_attainment(pops[pops != "fin"], sample_sizes, mode='world', gene=by_gene, rare=only_rare, hom=hom)
    save(world_attn_curve, file=args[4])
  } 
  else {
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

