elif args.gene_level:
	if args.output:
		with open(args.output, 'a') as out:
			out.write('GENE\tSYMBOL\tNUM_LoFs\tTEST\tP-VALUE\tCAC_Adj\t%s\tCAN_Adj\t%s\n' % ('\t'.join('CAC_%s' % pop for pop in pops), '\t'.join('CAN_%s' % pop for pop in pops)))
	num_biased_genes = 0
	gene_pop_ac = defaultdict(lambda: defaultdict(int)) # aggregate allele count for LoFs in a given gene
	gene_pop_an = defaultdict(lambda: defaultdict(int))
	gene_an = defaultdict(int)
	gene_num_lofs = defaultdict(int)
	gene_map = {}	# map gene to symbol
	while vcf.reading:
		vcf.read_line()
		genes = set()
		acs = map(float, vcf.info_field['AC_Adj'].split(','))
		an = float(vcf.info_field['AN_Adj'])
		alts = vcf_parser.alt_allele.split(',')
		lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)
		pop_acs = dict([(pop, map(float, vcf.info_field['AC_%s' % pop].split(','))) for pop in pops]) # maps pop to list of ACs
		pop_ans = dict([(pop, float(vcf.info_field['AN_%s' % pop])) for pop in pops]) # maps pop to float AN
				
		for allele_num in get_set_from_annotation(lof_annotations, 'ALLELE_NUM'):
			af = acs[int(allele_num) - 1]/an # global allele frequency for allele_num
			if af > AF_CUTOFF or af < LOW_AF_CUTOFF: continue
			snp = len(alts[int(allele_num) - 1]) == len(fields[header['REF']])
			this_alt_allele = filter_annotation(lof_annotations, 'ALLELE_NUM', allele_num) # list of annotations for allele_num
			# loop over genes associated with current alt allele
			for gene in get_set_from_annotation(this_alt_allele, 'Gene'):
				genes.add(gene)
				gene_num_lofs[gene] += 1
				gene_map[gene] = get_feature_from_annotation(filter_annotation(this_alt_allele, 'Gene', gene), 'SYMBOL')
				for pop in POPS:
					pop_ac = pop_acs[pop][int(allele_num) - 1]
					gene_pop_ac[gene][pop] += pop_ac
		
		for gene in genes:
			gene_an[gene] += an # global chromosome number for each gene aggregated over LoFs
			for pop in pops:
				gene_pop_an[gene][pop] += pop_ans[pop]

	# run chisquare/multinomial exact test for each gene as a whole
	for gene in gene_an:
		sig = False
		an_adj = gene_an[gene]
		prob_params = [gene_an[gene][pop]/an_adj for pop in pops]
		n = sum(gene_ac[gene][pop] for pop in pops)
		expected = map(lambda x:x*n, prob_params)
		observed = [gene_ac[gene][pop] for pop in pops]
		if sum(imap(lambda x: x<5, expected)) > 0:
			m = Multinomial(prob_params)
			p = m.exact_test(observed)
			exact = True
		else:
			chisq, p = chisquare(observed, f_exp=expected)
			exact = False

		if p < args.sig_thresh:
			sig = True
			num_biased_genes += 1
			if args.output:
				line = [gene, gene_map[gene], gene_num_lofs[gene]]
				line.append('MULTINOMIAL_EXACT' if exact else 'CHI_SQUARE')
				line.append(p)
				line.append(sum(gene_pop_ac[gene][pop] for pop in pops))
				line.extend(map(int, [gene_pop_ac[gene][pop] for pop in pops]))
				line.append(int(gene_an[gene]))
				line.append(sum(gene_pop_an[gene][pop] for pop in pops))
				with open(args.output, 'a') as out:
					out.write('%s\n' % '\t'.join(map(str,line)))
	print num_biased_genes


			