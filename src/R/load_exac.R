exac_all <- get(load('R_data/exac_pass.RData'))
exac_lof <- subset(exac, lof == 'HC' & is.na(lof_flags) & use)
pops <- c('afr', 'fin', 'nfe', 'amr', 'sas', 'eas')
pop_sizes <- apply(exac_lof[,paste("an",pops,sep="_")], 2, max)/2 # in people number
world_pop_sizes <- data.frame(
		row.names=c("eas", "sas", "nfe", "Middle Eastern", "afr", "amr", "Oceanic", "DiverseOther", "AfricanEuropeanAdmixed"),
		world=c(1932, 2085, 1145, 410,  1022, 529, 38, NA, NA)
)

