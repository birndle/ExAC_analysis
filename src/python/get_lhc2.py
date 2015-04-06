import pandas as pd
df = pd.io.parsers.read_table("exac_table.txt", index_col=0)
lof_subset = df.loc[df.use & df.lof_use,:] 

def count_lofs(args, lof_df):
	vcf_header = sp.check_output(['tabix', '-H', args.vcf]).strip().split('\n')
	header = VCF_Parser(vcf_header)
	peeps = header.sample_names
	header = header.header
	lof_counts = defaultdict(lambda: defaultdict(int)) # gene -> individual -> count
	hom_counts = defaultdict(lambda: defaultdict(int))
	haplotype_counts = defaultdict(lambda: defaultdict(int))

	for idx, row in lof_df.iterrows():
		region = '%s:%s-%s' % (row['chrom'], row['pos'], row['pos'])
		allele = row['alt']
		# genes = row['gene'] # check that table has this field
		vcf = sp.check_output(['tabix', args.full_vcf, region]).strip().split('\n')
		for line in vcf: # should actually be only one line
			line = line.split('\t')
			alt_alleles = line[header['ALT']].split(',')
			a_num = alt_alleles.index(allele) + 1
			for peep in peeps:
				gt = line[header[peep]]
				lofs = sum(1 for a in map(int, gt.split('/')) if a == a_num)
				for gene in genes:
					lof_counts[gene][peep] += lofs
				# how to count up haplotypes?





	# PSEUDOCODE
	# for line in sites.vcf
		# if LoF AND other cutoffs:
			# get position and overlapping gene(s)
			# tabix full.vcf position
			# aggregate # of LoF alleles per individual 
			# aggregate # of haplotypes per individual
				# THINK of a better way to count up HNs. Want it to approach 2 asymptotically

""" Return usable LoF alleles for a given loci in the sites vcf """
def use(vcf):
	lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)
	# APPLY AN and AF FILTERS HERE
	if lof_annotations:
		return True
	# lofs = filter_annotation_list(lof_annotations, 'Gene', [gene], True)
	# lof_alleles = get_set_from_annotation(lofs, 'ALLELE_NUM')

