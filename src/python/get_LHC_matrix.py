# CHECK OUTCOMES OF JOBS
# for chunk in ./*; do echo $(cat $chunk | head -n 18 | tail -n 1) $chunk; done
# ASSUMING ALL FAILURES ARE DUE TO MEMORY:
# for chunk in ./*; do echo $(cat $chunk | head -n 18 | tail -n 1) $chunk; done | grep TERM | awk ' BEGIN{FS="/"} {print $2}' > failed_chunks.txt

def split_region(r, num_chunks=20):
	from numpy import linspace
	c = r[0]
	segments = linspace(int(r[1])-1, int(r[2]), num_chunks, dtype=int)
	chunks = []
	for i in range(len(segments)-1):
		chunk = [c, segments[i]+1, segments[i+1]]
		chunks.append(chunk)
	return [map(str, r) for r in chunks]



def main(args):
	gene_positions = get_map(args.genes, strip=True, one_to_one=True)
	lof_counts = defaultdict(lambda: defaultdict(int))
	lhc = defaultdict(lambda: defaultdict(int))
	vcf_header = sp.check_output(['tabix', '-H', args.vcf]).strip().split('\n')
	indivs = VCF_Parser(vcf_header).sample_names
	hn = defaultdict(lambda: defaultdict(int))
	start = time.time()
	# get LoF counts for each individual for each gene
	for gene in gene_positions:
		region = gene_positions[gene]
		if args.split_genes:
			regions = split_region(region)
		else:
			regions = [region]
		regions = [map(str, r) for r in regions]

		for region in regions:
			print gene, '%s:%s-%s' % tuple(region)

			# get LoF positions for given gene
			vcf_sites_section = sp.check_output(['tabix', '-h', args.vcf_sites, '%s:%s-%s' % tuple(region)]) 
			vcf_sites_section = vcf_sites_section.strip().split('\n')
			lofs_pos, alts = get_LoFs_for_gene(vcf_sites_section, gene)
			if not lofs_pos: continue
			# get LoF counts for given gene for each individual
			vcf = sp.check_output(['tabix', '-h', args.vcf, '%s:%s-%s' % tuple(region)])
			vcf = VCF_Parser(vcf.strip().split('\n'))
			if vcf.line.startswith('#'): continue
			while vcf.reading:
				var = (vcf.get_field('CHROM'), vcf.get_field('POS'))
				if not var in lofs_pos: 
					vcf.read_line()
					continue
				allele_nums = alts[lofs_pos.index(var)]
				for indiv in vcf.sample_names:
					gt = vcf.get_genotype_for_individual(indiv)
					if gt == './.': continue
					hits = sum(map(lambda a:a in allele_nums, gt.split('/')))
					if not hn[gene][indiv]: hn[gene][indiv] += 2
					if hits == 2:
						lhc[gene][indiv] = hits
					lof_counts[gene][indiv] += hits
				vcf.read_line()

	# convert LoF counts to LHC
	for gene in lof_counts:
		for indiv in indivs:
			if lhc[gene][indiv] == 2: continue
			if lof_counts[gene][indiv] < 2:
				lhc[gene][indiv] = lof_counts[gene][indiv]
			else:
				lhc[gene][indiv] = 1 + (1-0.5**(lof_counts[gene][indiv]-1))

	if args.output:
		with open(args.output, 'a') as o:
			o.write('\t'.join(indivs)+'\n')
			for gene in lhc:
				line = ['%s,%s' % (lhc[gene][i], hn[gene][i]) for i in indivs]
				o.write('%s\t%s\n' % (gene, '\t'.join(line)))

# python /humgen/atgu1/fs03/birnbaum/ExAC_analysis/python_scripts/get_LHC_matrix.py --genes LHC/genes/chunk45 --vcf_sites /humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/ExAC.r0.3.sites.vep.vcf.gz --vcf /humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/exac_all.vcf.gz --output LHC/tables/chunk45

# pull LoF variants from sites vcf
def get_LoFs_for_gene(vcf, gene):
	vcf = VCF_Parser(vcf)
	if vcf.line.startswith('#'): return None, None
	positions = []
	alts = []
	while vcf.reading:
		if not vcf.annotations:
			vcf.read_line()
			continue
		# filter annotation list down to only HC_LoFs that don't carry the flags 'PHYLOCSF_WEAK' or 'PHYLOCSF_UNLIKELY_ORF'
		lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)
		lofs = filter_annotation_list(lof_annotations, 'Gene', [gene], True)
		if not lofs: 
			vcf.read_line()
			continue
		lof_alleles = get_set_from_annotation(lofs, 'ALLELE_NUM')
		chrom = vcf.get_field('CHROM')
		pos = vcf.get_field('POS')
		positions.append((chrom, pos))
		alts.append(lof_alleles)
		vcf.read_line()
	return positions, alts


if __name__ == '__main__':
	import argparse
	import os
	import sys
	sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/python_scripts')
	from utilities import get_map, parse_commandline_argument
	from vcf_parsing import VCF_Parser
	from parsing_exac_vcf import filter_annotation_list, get_set_from_annotation, filter_annotation 
	import subprocess as sp
	from collections import defaultdict
	from itertools import imap
	import time
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf', dest='vcf') 
	parser.add_argument('--vcf_sites', dest='vcf_sites')
	parser.add_argument('--genes', dest='genes')
	parser.add_argument('--parallelize', dest='parallelize', action='store_true')
	parser.add_argument('--gene_dir', dest='gene_dir')
	parser.add_argument('--chunk_genes', dest='chunk_genes', action='store_true')
	parser.add_argument('--log_dir', dest='log_dir')
	parser.add_argument('--output_dir', dest='output_dir')
	parser.add_argument('--output', dest='output')
	parser.add_argument('--chunk_list', dest='chunk_list')
	parser.add_argument('--chunk_size', dest='chunk_size', type=int)
	parser.add_argument('--split_genes', dest='split_genes', action='store_true')
	parser.add_argument('--merge_chunks', dest='merge_chunks', action='store_true')
	parser.add_argument('--unmerged', dest='unmerged', help='Directory of directories of LHC matrices that need to be merged together.')
	parser.add_argument('--extract_from', dest='extract_from')
	args = parser.parse_args()

	if args.chunk_genes:
		if args.chunk_list:
			genes = {}
			GENES = '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/LHC/sub_chunks/genes'
			with open(args.chunk_list, 'r') as chunks:
				for chunk in chunks:
					f = os.path.join(GENES, chunk.strip())
					with open(f, 'r') as i:
						for line in i:
							line = line.strip().split()
							genes[line[0]] = line[1:]
		else:
			genes = get_map(args.genes, strip=True, one_to_one=True)
		chunk = 1
		i = 0
		out = []
		for gene in genes:
			out.append('%s\t%s' % (gene, '\t'.join(genes[gene])))
			i += 1
			if i == args.chunk_size:
				with open(os.path.join(args.gene_dir, 'chunk%s' % chunk), 'a') as o:
					o.write('\n'.join(out))
				out = []
				i = 0
				chunk += 1
		if out:
			with open(os.path.join(args.gene_dir, 'chunk%s' % chunk), 'a') as o:
				o.write('\n'.join(out))
	
	elif args.parallelize:
		SCRIPT = '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/python_scripts/get_LHC_matrix.py'
		VCF_SITES = '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/ExAC.r0.3.sites.vep.vcf.gz'
		VCF_FULL = '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/exac_all.vcf.gz'
		if not (args.gene_dir and args.output_dir and args.log_dir):
			print >> sys.stderr, "Arguments \"--gene_dir\", \"--output_dir\", and \"--log_dir\" needed for paralleization. Exiting."
			sys.exit(1)

		# use args.chunk_list to subset args.gene_dir down to a list of chunks
		if args.chunk_list:
			with open(args.chunk_list, 'r') as i:
				chunks = [os.path.join(args.gene_dir, chunk.strip()) for chunk in i]
		else:
			chunks = parse_commandline_argument(args.gene_dir)
		for chunk in chunks:
			job_name = chunk.split('/')[-1]
			bsub_cmd = ['bsub', '-o', os.path.join(args.log_dir, job_name), '-q', 'hour', '-W', '1:00', '-J', job_name, '-R', 'rusage[mem=10]']
			output = os.path.join(args.output_dir, job_name)
			cmd = 'python %s --genes %s --vcf_sites %s --vcf %s --output %s' % (SCRIPT, chunk, VCF_SITES, VCF_FULL, output)
			if args.split_genes:
				cmd += ' --split_genes'
			bsub_cmd.append(cmd)
			# print ' '.join(bsub_cmd)
			sp.call(bsub_cmd)

	# merge multi-gene LHC tables that result from chunking up the gene list into smaller chunks of genes
	elif args.merge_chunks:
		if not (args.unmerged and args.output):
			print >> sys.stderr, "Arguments \"--unmerged\" and \"--output\" needed for merging. Exiting"
			sys.exit(1)
		merged = open(args.output, 'a')
		vcf_header = sp.check_output(['tabix', '-H', '/humgen/atgu1/fs03/birnbaum/ExAC_analysis/vcfs/exac_all.vcf.gz']).strip().split('\n')
		indivs = VCF_Parser(vcf_header).sample_names
		merged.write('\t'.join(indivs)+'\n')
		forks = [f for direc in parse_commandline_argument(args.unmerged) for f in direc]
		for fork in forks:
			skip_header = True
			with open(fork, 'r') as i:
				for line in i:
					if skip_header:
						skip_header = False
						continue
					merged.write(line)
		merged.close()

	elif args.extract_from:
		if not (args.output_dir):
			print >> sys.stderr, "Argument \"--output_dir\" needed for extracting. Exiting"
			sys.exit(1)
		header = False
		lhc_out = open(os.path.join(args.output_dir, 'LHC.table'), 'a')
		hn_out = open(os.path.join(args.output_dir, 'HN.table'), 'a')
		with open(args.extract_from, 'r') as i:
			for line in i:
				if not header:
					header = True
					lhc_out.write(line)
					hn_out.write(line)
					continue
				line = line.strip().split()
				lhc_line = [line[0]] + [x.split(',')[0] for x in line[1:]]
				hn_line = [line[0]] + [x.split(',')[1] for x in line[1:]]
				lhc_out.write('%s\n' % '\t'.join(lhc_line))
				hn_out.write('%s\n' % '\t'.join(hn_line))
		lhc_out.close()
		hn_out.close()

	else:
		main(args)

