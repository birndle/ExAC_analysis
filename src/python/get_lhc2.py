import os, sys
sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/src/python')

# genes = list of genes by which to subset the sites vcf 
def load_sites(genes):
	import pandas as pd
	SITES = "/humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/exac_sites_table.txt"
	print 'Loading in ExAC sites table ..'
	df = pd.io.parsers.read_table(SITES, index_col=0)
	lof_subset = df.loc[df.use & df.lof_use,:]
	print 'Done.'
	if genes:
		pattern = '|'.join(genes.split(','))
		lof_subset = lof_subset.loc[lof_subset.gene.str.contains(pattern)]
	return(lof_subset)


def get_genes():
	path = '/humgen/atgu1/fs03/birnbaum/gencode_annotations/gencode_all_genes.txt'
	with open(path, 'r') as f:
		genes = [line.strip() for line in f]
	return(genes)


def count_lofs(args, lof_df):
	import subprocess as sp
	from vcf_parsing import VCF_Parser
	from collections import defaultdict

	vcf_header = sp.check_output(['tabix', '-H', args.vcf]).strip().split('\n')
	vcf_header = VCF_Parser(vcf_header)
	peeps = vcf_header.sample_names
	header = vcf_header.header
	lof_counts = defaultdict(lambda: defaultdict(int)) # gene -> individual -> count
	is_hom = defaultdict(lambda: defaultdict(int))
	haplotype_counts = defaultdict(lambda: defaultdict(int))

	gt_idx = None
	print "Counting LoFs .."
	for idx, row in lof_df.iterrows():
		region = '%s:%s-%s' % (row['chrom'], row['pos'], row['pos'])
		print region
		allele = row['alt']
		genes = row['gene'].split(',') # check that table has this field
		vcf = sp.check_output(['tabix', args.vcf, region]).strip().split('\n')
		for line in vcf: # should actually be only one line
			line = line.split('\t')
			if not gt_idx:
				gt_idx = line[header['FORMAT']].split(':').index('GT')

			if line[header['POS']] != row['pos']: 
				continue
				# this happens when an indel overlaps the current position and thus is returned by tabix,
				# even though it is not at the position we were are currently scanning from the sites table

			alt_alleles = line[header['ALT']].split(',')
			a_num = alt_alleles.index(allele) + 1
			for peep in peeps:
				gt = line[header[peep]].split(':')[gt_idx]
				if gt == './.':
					continue

				lofs = sum(1 for a in map(int, gt.split('/')) if a == a_num)
				for gene in genes:
					lof_counts[gene][peep] += lofs
					if lofs == 2:
						is_hom[gene][peep] = 1
					haplotype_cts[gene][peep] += 2

	genes = get_genes()
	def write_results(d, which):		
		# write counts to output file
		output = os.path.join(args.output_dir, which, '%s.txt' % args.idx)
		o = open(output, 'a')
		o.write('\t'.join(peeps) + '\n') # column names
		lines = ['%s\t' % gene + '\t'.join(map(str, [d[gene][peep] for peep in peeps])) for gene in genes]
		o.write('\n'.join(lines))
		o.close()

	print 'Finished counting. Now writing results.'
	print 'Writing LoF counts table ..'
	write_results(lof_counts, 'lof_cts')
	print 'Writing homozygote truth table ..'
	write_results(is_hom, 'has_hom')
	print 'Writing haplotype counts table ..'
	write_results(haplotype_counts, 'haplotype_cts')
	print 'Done.'


def main(args):
	import numpy as np
	lof_sites = load_sites(args.genes)
	print 'Subsetted table contains %s variant loci.' % lof_sites.shape[0]
	nrow = lof_sites.shape[0]
	chunks = np.linspace(0, nrow, args.chunks, dtype=int)
	subset = lof_sites.iloc[chunks[args.idx]:chunks[args.idx+1]]
	print 'This chunk contains %s variant loci.' % subset.shape[0]
	count_lofs(args, subset)

# e.g. in ExAC_analysis directory
# python src/python/get_lhc2.py --parallelize --output_dir lhc_test --log_dir lhc_test/logs --chunks <num_chunks> --test
def parallelize(args):
	import subprocess as sp

	SCRIPT = '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/src/python/get_lhc2.py'
	VCF_FULL = '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/exac_all.vcf.gz'
	if not (args.output_dir and args.log_dir and args.chunks):
		print >> sys.stderr, "--output_dir, --log_dir, & --chunks parameters needed for paralleization. Exiting."
		sys.exit(1)

	for i in range(0, args.chunks-1): # [0, args.chunks-1)
		job = 'lhf%s' % i
		bsub_cmd = ['bsub', '-oo', os.path.join(args.log_dir, job), '-q', 'hour', '-W', '4:0w0', '-J', job, '-R', 'rusage[mem=20]']
		cmd = 'python %s --vcf %s --chunks %s --idx %s --output_dir %s' % (SCRIPT, VCF_FULL, args.chunks, i, args.output_dir)
		bsub_cmd.append(cmd)
		sp.call(bsub_cmd)

		if args.test:
			break


# aggregate the count tables coming from different chunks of the sites vcf
def merge(args):
	from utilities import listdir_fullpath
	import pandas as pd
	chunk_dir = os.path.join(args.input_dir, args.mode)
	df = None
	for chunk in listdir_fullpath(chunk_dir):
		if not df:
			df = pd.io.parsers.read_table(chunk)
			continue
		
		new_df = pd.io.parsers.read_table(chunk)
		if args.mode == 'lof_cts' or args.mode == 'haplotype_cts':
			df = df + new_df
		elif args.mode == 'has_hom':
			df = ((df == 1) or (new_df == 1))*1 
	
	output = os.path.join(args.input_dir, '%s.txt' % args.mode)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf', help='Full ExAC vcf, with genotypes.', dest='vcf')
	parser.add_argument('--parallelize', action='store_true', dest='parallelize')
	parser.add_argument('--log_dir', dest='log_dir', help='Path to log file for bsub job.')
	parser.add_argument('--chunks', dest='chunks', help='Integer, number of chunks to split sites df into.', type=int)
	parser.add_argument('--idx', dest='idx', help='Which chunk to process. Gets passed for each chunk during parallelization.', type=int)
	parser.add_argument('--output_dir', dest='output_dir', help='Where to write processed chunks.')
	parser.add_argument('--input_dir', dest='input_dir', help='Where chunks have been written.')
	parser.add_argument('--genes', dest='genes', help='Genes by which to subset', required=False)
	parser.add_argument('--merge', dest='merge', action='store_true')
	parser.add_argument('--mode', dest='mode', help='\'lof_cts\', \'has_hom\', or \'haplotype_cts\'. Only use with --merge.')
	parser.add_argument('--test', dest='test', action='store_true')
	args = parser.parse_args()
	
	if args.parallelize:
		parallelize(args)
	elif args.merge:
		merge(args)
	else:
		main(args)


