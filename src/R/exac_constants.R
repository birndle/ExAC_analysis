# ExAC constants file
# Use by source('exac_constants.R')

options(stringsAsFactors=FALSE)
library(plyr)

# The colors!
color_amr = k_amr = '#32C832'
color_eur = k_eur = '#3299CC'
color_afr = k_afr = '#CD2626'
color_sas = k_sas = '#643200'
color_eas = k_eas = '#FF9B00'
color_oth = k_oth = '#CDC9C9'

color_syn = k_syn = '#AAAAAA'
color_mis = k_mis = '#FF6103'
color_lof = k_lof = '#9D1309'

# example usage: alpha(k_lof,.5) gives you a 50% transparent LoF maroon
alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

pops = c('AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS')

## Population sizes
# Populations and percent ancestries
# World in millions, all others actual
# From Laramie
populationAncestries <- data.frame(
  row.names=c("East Asian", "South Asian", "European", "Middle Eastern", "African", "Latino", "Oceanic", "DiverseOther", "AfricanEuropeanAdmixed"),
  world=c(1932, 2085, 1145, 410,  1022, 529, 38, NA, NA),
  kgenomes=c(523, 494, 514, 0, 691, 355, 0, NA, NA),
  esp=c(0, 0, 4298, 0, 2217,0,0, NA, NA),
  exac=c(4327, 8256, 36677, 0, 5203, 5789, 0,454, NA),  # Numbers from Monkol on 3/4/15
  pcgcGwas=c(5219, 0, 116766, 0, 0, 0, 0, NA, NA),   # PGC Published 4 (SCZ, BIP, MDD, ADHD)
  ptsd=c(106, 0, 8393, 0, NA, 829, 0, 1295, 9845)
)

# Loading data
data_url = 'http://www.broadinstitute.org/~konradk/exac_data/'
open_file_exac = function(fname) {
  con = gzcon(url(paste0(data_url, fname)))
  dat = readLines(con)
  return(read.delim(textConnection(dat), header=T))
}
open_file_exac_download = function(fname) {
  download.file(paste0(data_url, fname), fname)
  return(read.delim(fname, header=T))
}

load_exac_data = function(type='', hgmd=TRUE, local=FALSE) { # This takes ~10 minutes (download, then ~5 mins of processing)
  transcripts = ''
  print('Loading data...')
  if (type == 'canonical') {
    transcripts = '.canonical'
  }
  if (hgmd) {
    data_format = '.full'
  } else {
    data_format = '.popmax.clinvar'
  }
  fname = paste0('ExAC_HC.0.3.final.vep', data_format, transcripts, '.table.gz')
  if (local) {
    exac = read.delim(fname, header=T)
  } else {
    exac = open_file_exac_download(fname)
  }
  print('Done! Processing data... Computing af, maf, indel, bases_inserted...')
  colnames(exac) = tolower(colnames(exac))
  exac$af = 0.0
  exac$af[exac$an_adj > 0] = as.numeric(exac$ac_adj[exac$an_adj > 0]) / as.numeric(exac$an_adj[exac$an_adj > 0])
  exac$maf = pmin(exac$af,1-exac$af) # minor allele frequency
  exac$indel = nchar(exac$ref) != nchar(exac$alt)
  exac$bases_inserted = nchar(exac$alt) - nchar(exac$ref)
  exac$transition = (exac$ref == 'A' & exac$alt == 'G') | (exac$ref == 'G' & exac$alt == 'A') | (exac$ref == 'C' & exac$alt == 'T') | (exac$ref == 'T' & exac$alt == 'C')
  exac$transition[exac$indel] = NA
  
  print('Done! Calculating bad regions of the genome and what to use...')
  resolution = 1000
  intervals=(0:(250000000/resolution))*resolution
  exac$pos_bin = cut(exac$pos, intervals)
  exac$bin_start = as.numeric( sub("\\((.+),.*", "\\1", exac$pos_bin))
  exac$bin_end = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", exac$pos_bin))
  allelic_state = count(subset(exac, select=c(chrom, pos)))
  multiallelics = subset(allelic_state, freq > 3, select=c(chrom, pos))
  multiallelics$pos_bin = cut(multiallelics$pos, intervals)
  multiallelics$bin_start = as.numeric( sub("\\((.+),.*", "\\1", multiallelics$pos_bin))
  multiallelics$bin_end = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", multiallelics$pos_bin))
  multiallelic_counts = count(multiallelics, vars = c('chrom', 'bin_start', 'bin_end'))
  bad_sectors = subset(head(multiallelic_counts[order(multiallelic_counts$freq, decreasing = T),], 10), select=c(chrom, bin_start))
  bad_sectors$bad = TRUE
  exac = merge(exac, bad_sectors, by=c('chrom', 'bin_start'), all.x=T)
  exac$use = exac$an_adj > .95*max(exac$an_adj, na.rm=TRUE) & exac$ac_adj > 0 & exac$filter=='PASS' & is.na(exac$bad)
  exac$lof_use = !is.na(exac$lof) & exac$lof == 'HC' & is.na(exac$lof_flags)
  
  print('Parsing SIFT and PolyPhen scores...')
  # for the idea behind these named group regexes, see http://stackoverflow.com/a/2969666/3806692
  sift_regex = "([a-z]*)\\(([0-9\\.]*)\\)"
  exac$sift_word = sub(sift_regex, "\\1", exac$sift)
  exac$sift_score = as.numeric(sub(sift_regex, "\\2", exac$sift)) 
  polyphen_regex = "([a-z_]*)\\(([0-9\\.]*)\\)"
  exac$polyphen_word = sub(polyphen_regex, "\\1", exac$polyphen) # parse words
  exac$polyphen_score = as.numeric(sub(polyphen_regex, "\\2", exac$polyphen)) # parse scores
  
  print('Separating functional categories...')
  exac$category = exac$consequence
  exac$category[exac$consequence=='inframe_insertion'] = 'inframe_indel'
  exac$category[exac$consequence=='inframe_deletion'] = 'inframe_indel'
  exac$category[exac$category=="\\N"] = NA

  print('Computing allele frequencies...')
  exac$af_popmax = exac$ac_popmax / exac$an_popmax
  exac$af_popmax[exac$an_popmax == 0 | is.na(exac$an_popmax) | is.na(exac$ac_popmax)] = 0.0
  exac$af_global = exac$ac_adj / exac$an_adj
  exac$af_global[exac$an_adj == 0 | is.na(exac$an_adj) | is.na(exac$ac_adj)] = 0.0
  exac$singleton = exac$ac_adj == 1
  
  print('Done! Here you go.')
  return(exac)
}

load_constraint_data = function() {
  constraint = open_file_exac('fordist_cleaned_exac_r03_march16_z_data.txt.gz')
  constraint$feature = sapply(constraint$transcript, function(x) { strsplit(x, '.', fixed=TRUE)[[1]][[1]] })
  
  #boundaries = c(0, 0.25, 0.75, 1)
  boundaries = c(0, 0.25, 0.5, 0.75, 1)
  constraint$lof_cut = cut(constraint$corrected_lof, quantile(constraint$corrected_lof, probs=boundaries))
  constraint$mis_cut = cut(constraint$corrected_mis, quantile(constraint$corrected_mis, probs=boundaries))
  constraint$syn_cut = cut(constraint$corrected_syn, quantile(constraint$corrected_syn, probs=boundaries))
  
  # Confirmed to Kaitlin's criteria (4982 transcripts)
  constraint$lof_constrained = constraint$adj_exp_lof >= 10 & constraint$corrected_lof > 3.09 & constraint$corrected_syn < 3.09
  constraint$mis_constrained = constraint$corrected_mis > 3.09 & constraint$corrected_syn < 3.09
  return(constraint)
}

caf_pops = c('AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS')
pop_cafs = c('AFR_CAF', 'AMR_CAF', 'EAS_CAF', 'FIN_CAF', 'NFE_CAF', 'SAS_CAF')
# Ordering by code alphabetical with oth at the end
pop_colors = c(color_afr, color_amr, color_eas, alpha(color_eur, 0.5), color_eur, color_sas)

load_caf_data = function(type='') {
  if (type == 'transcript') {
    print('Loading transcript data...')
    caf = open_file_exac('ExAC_HC.0.3.final.vep.transcript.caf.gz')
  } else {
    print('Loading gene data...')
    caf = open_file_exac('ExAC_HC.0.3.final.vep.gene.caf.gz')    
  }
  print('Done!')
  
  maxes = t(apply(caf, 1, function(x) { 
    c(max(as.numeric(x[pop_cafs])),
      caf_pops[which.max(as.numeric(x[pop_cafs]))],
      mean(as.numeric(x[pop_cafs])), 
      median(as.numeric(x[pop_cafs])),
      min(as.numeric(x[pop_cafs]))
      )
  }))
  
  caf$max_caf = as.numeric(maxes[,1])
  caf$max_caf_pop = maxes[,2]
  caf$average_caf = as.numeric(maxes[,3])
  caf$median_caf = as.numeric(maxes[,4])
  caf$min_caf = as.numeric(maxes[,5])
  caf$delta_max = caf$max_caf/caf$average_caf
  return(caf)
}

load_gene_lists = function() {
#   list_metadata = read.table(textConnection("
# filename|display|xval
# universe|All genes|-2
# all_ar|Autosomal recessive|-3
# all_ad|Autosomal dominant|-4
# core_essentials_hart|Essential in culture|-5
# haploinsufficient|Haploinsufficient|-6
# mgi_essential|Essential in mice|-7
# fda_approved_drug_targets|Drug targets|-8
# grep_or|Olfactory receptors|-9
# gwascatalog|Near GWAS hits|-10
# homozygous_lof_tolerant_twohit|Homozygous LoF tolerant|-11
# "),sep='|',header=TRUE)
  
  if (!file.exists('gene_lists')) {
    system('git clone git@github.com:macarthur-lab/gene_lists.git')
  } else {
    print('Using locally available gene_lists... Run git pull to update.')
  }
  gene_lists = ldply(Sys.glob('gene_lists/lists/*'), function(x) {
#     gene_list = read.table(paste('gene_lists/lists/', x, '.tsv',sep=''))
    gene_list = read.table(x)
    names(gene_list) = 'symbol'
    gene_list$list = strsplit(strsplit(x, '/', fixed=T)[[1]][3], '.', fixed=T)[[1]][1]
    gene_list
  })
  
  all_genes = dcast(gene_lists, symbol ~ list, fun.aggregate = function(x) { length(x) >= 1 })
}
