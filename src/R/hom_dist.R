hom_dist <- function(pop, hom_df) {
  pop_homs <- data[,paste("hom",pop,sep="_")]
  pop_sizes <- round(data[,paste("an",pop,sep="_")]/2)
  hom_prev <- pop_homs/pop_sizes
  breaks <- c(-Inf, 0.0005, 0.001, 0.005, 0.01, 0.05, Inf)
  num_bins <- length(breaks)-1
  bins <- cut(pop_af, breaks, labels=1:num_bins)
}