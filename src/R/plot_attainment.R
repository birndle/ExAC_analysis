load_samples <- function(path) {
  samples <- list.files(path, full.names=TRUE)
  prof <- get(load(samples[1]))
  if (length(samples) == 1) {
    return(prof)
  }
  for (i in 2:length(samples)) {
    prof <- prof + get(load(samples[i]))   
  }
  return(prof/length(samples))
}

#  conds is a string, e.g. "gene_all" or "lof_hom"
draw_plot <- function(conds, title, y, ylim, xlim=130000/2, diverse=F, legloc="topleft") {
  source("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/Rscripts/exac_colors.R")
  prefix <- sprintf("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/attainment/curves/pops/%s", conds)
  paths <- sapply(pops, function(pop) { sprintf("%s/%s", prefix, pop) })
  pop_curves <- sapply(paths, load_samples)
  pop_curves[pop_curves==0] <- -Inf
  xax=seq(from=1000,to=121000,by=500)
  pop_curves <- pop_curves[1:length(xax),]
  
  if (diverse) {
    unif <- load_samples(sprintf("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/attainment/curves/all/uniform/%s",conds))
    unif[unif==0] <- -Inf
    unif <- unif[1:length(xax)]
    exac_proportioned <- load_samples(sprintf("/Users/birnbaum88/Desktop/Macarthur/ExAC_analysis/attainment/curves/all/exac/%s",conds))
    exac_proportioned[exac_proportioned==0] <- -Inf
    exac_proportioned <- exac_proportioned[1:length(xax)]
    pop_curves <- cbind(pop_curves, unif, exac_proportioned)
    pop_colors <- c(pop_colors, "purple", "pink")
    pops <- c(pops, "Uniformly portioned", "ExAC portioned")
  }
  
  matplot(xax/2, pop_curves, col=pop_colors, type='l', lty=1, lwd=3, xlim=c(0,xlim), ylim=c(0,ylim), main=title,xlab="Sample size in people",ylab=y)
  legend(legloc, pops, fill=pop_colors, cex=0.6, ncol=1)
}


