pops <- c('afr', 'fin', 'nfe', 'amr', 'sas', 'eas')
colors <- list()
colors['amr'] = '#32C832'
colors['nfe'] = '#3299CC'
colors['afr'] = '#CD2626'
colors['sas'] = '#643200'
colors['eas'] = '#FF9B00'
colors['oth'] = '#CDC9C9'
colors['fin'] = '#003580'
pop_colors <- sapply(pops, function(pop) {return(colors[[pop]])})
