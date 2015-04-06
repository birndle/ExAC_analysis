from subprocess import Popen, PIPE, STDOUT
import time
from collections import defaultdict
import os

"""
Input: coverageBed output file name OR list of lines in coverage files, extracted by "get_profile_for_region" 
"""
def parse_pileup(pileup, d=True, sort=None, strand_dict=None, bed_dict=None, test=False):	
	if type(pileup) == str:
		# read pileup file 
		with open(pileup, 'r') as p:
			lines = [line.strip().split('\t') for line in p]
	else:
		lines = [line.split('\t') for line in pileup if not line == '']

	if sort:
		lines = sorted(lines, key=lambda x: float(x[1])/float(x[2]))

	if test:
		print [l[:3] for l in lines] == test
		return

	if d:
		gene = None
		genes = []
		genome = defaultdict(dict)
		strands = []
		feats = []
		beds = []

		for line in lines:
			bed = tuple(line[:3])
			chrom = line[0]
			lb = int(line[1])
			pos = int(line[3])
			cvg = int(line[4])
			genome[chrom][lb+pos] = cvg
			if pos == 1:
				# retrieve information about the current genomic interval .. strand? feature?
				if bed_dict:
					feat = bed_dict[bed]
					if strand_dict:
						strand = strand_dict[feat]
						strands.append(strand)
					feats.append(feat)
					beds.append(bed)
				if gene:
					# aggregate dictionaries in a list
					genes.append(gene)
				gene = {}
			# populate a dictionary for each region of coverage, mapping pos to cvg
			gene[pos] = cvg
		# add last exon
		genes.append(gene)
		return genes, genome, strands, feats, beds
	else:
		# return read counts for each interval in pileup
		return [int(line[3]) for line in lines]

"""
Utility for running quicking coverage jobs.
Accepts either a list of bed tuples (e.g. ('X', 123435, 643155)), or the name of a bedfile
"""

def get_profile_for_region(regions, bam, d=True, output=None):
	if output and os.path.isfile(output):
		print '%s already exists. Not computing coverage.' % output
		return

	sample = bam.split('/')[-1].split('_')[0]
	# If given python list of bed tuples
	if type(regions) == list:
		regions = ['\t'.join([str(x) for x in region]) for region in regions]
		regions = '\n'.join(regions) + '\n'

		# Call coverageBed via subprocess module, piping in bedfile-like-region-object
		args = ['coverageBed', '-d', '-split', '-abam', bam, '-b', 'stdin']
		if not d:
			del args[1]
		proc = Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
		# form_regions = 'chr%s:%s-%s' % tuple(regions.strip().split('\t'))
		# print 'Getting coverage for %s over %s ...' % (bam.split('/')[-1], form_regions)
		start = time.time()
		pileup = proc.communicate(regions)[0].split('\n')

	# If given a bed file
	elif type(regions) == str:
		args = ['coverageBed', '-d', '-split', '-abam', bam, '-b', regions]
		if not d:
			del args[1]
		proc = Popen(args, stdout=PIPE, stderr=STDOUT)
		print 'Getting coverage for %s over %s ...' % (bam.split('/')[-1], regions)
		start = time.time()
		pileup = proc.communicate()[0].split('\n')

	if output:
		with open(output, 'a') as o:
			for line in pileup:
				o.write('%s\n' % line)
		print 'Took %s minutes.' % ((time.time()-start)/60.0)
	else:
		print 'Took %s minutes.' % ((time.time()-start)/60.0)
		return pileup


def get_list_of_all_bams(which_bams='all'):
	if which_bams == 'all':
		bam_names = '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools/bams.txt'
	elif which_bams == 'filtered':
		bam_names = '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools/filtered_bams.txt'
	with open(bam_names, 'r') as f:
		bams = [line.strip() for line in f]
	return bams


def listdir_fullpath(direc):
	import os
	direcf = [os.path.join(direc, f) for f in os.listdir(direc)]
	for f in direcf:
		try:
			yield os.readlink(f)
		except OSError:
			yield f


def get_depth_dict(info='/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/the_data/depth_distributions/depth_by_sample.txt'):
	with open(info, 'r') as h:
		depth_dict = dict([(line.split()[0], int(line.split()[1])) for line in h])
	return depth_dict


def bed_to_tuple(bed):
	with open(bed, 'r') as f:
		beds = [line.strip().split() for line in f]
	beds = [(x[0], int(x[1]), int(x[2])) for x in beds]
	return beds


def parse_commandline_argument(arg, return_list=False):
	# if user inputs filename
	if os.path.isfile(arg):
		if return_list:
			return [arg]
		else:
			return arg
	# if user inputs directory
	elif os.path.isdir(arg):
		return [parse_commandline_argument(f) for f in listdir_fullpath(arg)]
	elif len(arg.split(',')) > 1:
		return arg.split(',')
	elif not arg:
		return []

def parse_single_column_file(fi, typ='string'):
	with open(fi, 'r') as f:
		if not typ == 'string':
			if typ == 'int':
				lines = [int(line.strip()) for line in f]
			elif typ == 'float':
				lines = [float(line.strip()) for line in f]
		else:
			lines = [line.strip() for line in f]
	return lines


def get_map(in_file, reverse=False, strip=False, one_to_one=True, return_list=False, header=False):
	if one_to_one:
		the_map = {}
	else:
		the_map = defaultdict(list)
	the_header=None
	with open(in_file, 'r') as f:
		for line in f:
			line = line.strip().split()
			if header:
				the_header = line
				header = False
				continue
			transcript = line[0]
			if strip:
				transcript = transcript.split('.')[0]
			info = line[1:]
			if len(info) == 1 and not return_list:
				info = info[0]
			if reverse:
				the_map[tuple(info)] = transcript
			elif one_to_one:
				the_map[transcript] = info
			else:
				the_map[transcript].append(info)
	if the_header:
		return the_map, the_header
	else:
		return the_map

def get_bamfile_for_sample(sample, filtered=False):
	if filtered:
		return '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/reformed_bams/filtered/%s.filtered.bam' % sample
	else:
		return '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/bamfiles/%s_uniquely_mapped_reads.bam.sorted.bam' % sample


def parse_read_count_table(table, strip=False):
	counts = defaultdict(dict)
	with open(table, 'r') as f:
		header = None
		for line in f:
			if not header:
				header = line.strip().split('\t')
				continue
			else:
				line = line.split('\t')
				feature = line[0]
				if strip:
					feature = feature.split('.')[0]
				for i in range(len(line[1:])):
					sample = header[i]
					counts[sample][feature] = int(line[i+1])
	return counts


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--bed', dest="bed", help='Bedfile containing regions to compute coverage over', required=False)
	parser.add_argument('--counts', dest="counts", help='Output file containing table of read counts', required=False)
	parser.add_argument('--features', dest="features", help='Input file containing names of features in bedfile', required=False)
	parser.add_argument('--bam', dest="bam", help='Input file containing list of bams. If want all bams, enter "all"', required=True)
	args = parser.parse_args()

	if args.bam == 'all':
		bams = get_list_of_all_bams()
	elif args.bam:
		bams = get_list_of_all_bams(args.bam)

	if args.features:
		features = get_list_of_all_bams(args.features)

	samples = [x.split('/')[-1].split('_')[0] for x in bams]	
	bed = bed_to_tuple(args.bed)
	# print samples
	results = []
	counts = parallelize_coverageBed(bed, bams)
	# print counts
	write_read_count_table(args.counts, counts, samples, features)













