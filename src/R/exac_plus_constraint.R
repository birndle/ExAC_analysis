# ExAC with constraint

source('exac_constants.R')
exac = load_exac_data('canonical')
constraint = load_constraint_data()

exac_constraint = merge(exac, constraint, by='Feature')
columns = c('CHROM', 'POS', 'REF', 'ALT', 'PolyPhen', 'SIFT', 'Feature', 'Gene', 
            'corrected_syn', 'corrected_mis', 'corrected_lof')
out_con = gzfile('exac_plus_constraint.txt.gz', 'w')
write.table(subset(exac_constraint, use & Consequence == 'missense_variant', select=columns), file=out_con, quote=F, row.names=F, sep='\t')
close(out_con)