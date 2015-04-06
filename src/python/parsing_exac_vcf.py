def main(args):

	if args.filter:
		vcf = VCF_Parser(args.vcf, subset=args.vcf_output)
		while vcf.reading:
			if not vcf.annotations:
				vcf.read_line()
				continue
			for a in vcf.annotations:
				if a['LoF'] == 'HC':
					vcf.write_line() # will only write if args.output != None
					break
			vcf.read_line()

	elif args.pop_spec:
		pops = ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
		# pops = ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS']
		vcf = VCF_Parser(args.vcf, subset=args.vcf_output)
		# DON'T aggregate across genes
		if args.variant_level:
			if args.pop_max:
				pop_max_vs_af = []
				while vcf.reading:
					an = float(vcf.info_field['AN_Adj'])
					if an == 0: 
						vcf.read_line()
						continue
					acs = map(float, vcf.info_field['AC_Adj'].split(','))
					pop_acs = dict([(pop, map(float, vcf.info_field['AC_%s' % pop].split(','))) for pop in pops])
					pop_ans = dict([(pop, float(vcf.info_field['AN_%s' % pop])) for pop in pops])
					lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)
					for allele_num in get_set_from_annotation(lof_annotations, 'ALLELE_NUM'):
						idx = int(allele_num)-1
						global_af = acs[idx]/an
						if global_af == 0: continue
						pop_max = 0
						max_pop = ''
						for pop in pops:
							pop_ac = pop_acs[pop][idx]
							pop_an = pop_ans[pop]
							if pop_an == 0: continue
							af = pop_ac/pop_an
							if af > pop_max:
								pop_max = af
								max_pop = pop
						pop_max_vs_af.append((global_af, pop_max, max_pop, acs[idx]))
					vcf.read_line()

				if args.output:
					with open(args.output, 'a') as out:
						for point in pop_max_vs_af:
							out.write('%s\n' % '\t'.join(map(str, point)))

			# look for variants that have AC distribution significantly differing from expectation
			elif args.chisquare:
				num_tests = 0
				candidates = []
				if args.output:
					with open(args.output, 'a') as out:
						out.write('CHROM\tPOS\tREF\tALT\tTEST\tP-VALUE\tAC_Adj\t%s\tAN_Adj\t%s\n' % ('\t'.join('AC_%s' % pop for pop in pops), '\t'.join('AN_%s' % pop for pop in pops)))
				while vcf.reading:
					sig = False
					ac_totals = vcf.info_field['AC_Adj'].split(',')
					an_adj = float(vcf.info_field['AN_Adj'])
					pop_acs = [vcf.info_field['AC_%s' % pop].split(',') for pop in pops] # CHANGE TO DICT
					pop_ans = dict([(pop, float(vcf.info_field['AN_%s' % pop])) for pop in pops])
					prob_params = np.array([pop_ans[pop]/an_adj for pop in pops])
					lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)

					for allele_num in get_set_from_annotation(lof_annotations, 'ALLELE_NUM'):
						ac_total = int(ac_totals[int(allele_num)-1])
						if ac_total == 0:
							continue
						observed_counts = [int(o[int(allele_num)-1]) for o in pop_acs]
						expected_counts = ac_total*prob_params
						# rule of thumb: check that no expected counts are less than 5 (source: Wikipedia) 
						if sum(imap(lambda x: x<5, expected_counts)) > 0:
							m = Multinomial(prob_params)
							p = m.exact_test(observed_counts)
							exact = True
						# otherwise, if sample size is too small, use an exact multinomial test
						else:
							chisq, p = chisquare(observed_counts, f_exp=expected_counts)
							exact = False

						if args.output and p < args.sig_thresh:
							line = map(vcf.get_field, ['CHROM', 'POS', 'REF'])
							line.append(vcf.get_alt(allele_num))
							line.append('MULTINOMIAL_EXACT' if exact else 'CHI_SQUARE')
							line.append(p)
							line.append(int(ac_total))
							line.extend(map(int, observed_counts))
							line.append(int(an_adj))
							line.extend(map(int, [pop_ans[pop] for pop in pops]))
							with open(args.output, 'a') as out:
								out.write('%s\n' % '\t'.join(map(str,line)))
							vcf.write_line()
						num_tests += 1
					vcf.read_line()

			elif args.loo:
				header = None
				p_values = defaultdict(list)
				ns = []
				lines = []
				while vcf.reading:
					# print vcf.line
					an = float(vcf.info_field['AN_Adj'])
					ref = vcf.ref_allele
					chrom = vcf.get_field('CHROM')
					pos = vcf.get_field('POS')
					acs = map(float, vcf.info_field['AC_Adj'].split(','))
					pop_acs = dict([(pop, map(float, vcf.info_field['AC_%s' % pop].split(','))) for pop in pops])
					pop_ans = dict([(pop, float(vcf.info_field['AN_%s' % pop])) for pop in pops])
					if an == 0 or not vcf.annotations: 
						vcf.read_line()
						continue
					lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)
					for allele_num in get_set_from_annotation(lof_annotations, 'ALLELE_NUM'):
						idx = int(allele_num)-1
						n = acs[idx]
						alt = vcf.alt_allele.split(',')[idx]
						i = 0
						for pop in pops:
							if pop_ans[pop] == 0 or sum(pop_ans[p] for p in pops if p != pop) == 0:
								p_val = 'N/A'
							elif args.lrt: 
								p_val, stat = likelihood_ratio_test([pop_acs[p][idx] for p in pops], [pop_ans[p] for p in pops], i, only_enriched=args.enriched)
							elif args.fisher:
								p_val = fishers_exact([pop_acs[p][idx] for p in pops], [pop_ans[p] for p in pops], i)
							p_values[pop].append(p_val)
							i += 1
						ns.append(n)
						lines.append([(chrom, pos, ref, alt), [pop_acs[pop][idx] for pop in pops], [pop_ans[pop] for pop in pops]])
					vcf.read_line()	

				if args.output:
					with open(args.output, 'a') as out:
						out.write('CHROM\tPOS\tREF\tALT\tAC_Adj\t%s\tAN_Adj\t%s\t%s\n' % ('\t'.join('AC_%s' % pop for pop in pops), '\t'.join('AN_%s' % pop for pop in pops), '\t'.join('P_%s' % pop for pop in pops)))
						for i in range(len(ns)):
							n = ns[i]
							line = lines[i]
							ps =  map(str, [p_values[pop][i] for pop in pops])
							line = ['\t'.join(map(str, line[0])), n, '\t'.join(map(str, line[1])), sum(line[2]), '\t'.join(map(str, line[2])), '\t'.join(ps)]
							out.write('%s\n' % '\t'.join(map(str, line)))

			elif args.pop:
				total_pop_an = 0
				total_background_an = 0
				total_pop_ac = 0
				total_background_ac = 0
				popAFvsglobalAF = []
				while vcf.reading:
					an = float(vcf.info_field['AN_Adj'])
					acs = map(float, vcf.info_field['AC_Adj'].split(','))
					pop_acs = dict([(pop, map(float, vcf.info_field['AC_%s' % pop].split(','))) for pop in pops])
					pop_ans = dict([(pop, float(vcf.info_field['AN_%s' % pop])) for pop in pops])
					pop_an = pop_ans[args.pop]
					background_an = sum(pop_ans[pop] for pop in pops if pop != args.pop)
					if pop_an == 0 or background_an == 0 or not vcf.annotations: 
						vcf.read_line()
						continue
					total_pop_an += pop_an
					total_background_an += background_an
					lof_annotations = filter_annotation_list(filter_annotation(vcf.annotations, 'LoF'), 'LoF_flags', ['PHYLOCSF_WEAK', 'PHYLOCSF_UNLIKELY_ORF'], False)
					for allele_num in get_set_from_annotation(lof_annotations, 'ALLELE_NUM'):
						idx = int(allele_num)-1
						pop_ac = pop_acs[args.pop][idx]
						total_pop_ac += pop_ac
						background_ac = sum(pop_acs[pop][idx] for pop in pops if pop != args.pop)
						total_background_ac += background_ac
						pop_af = pop_ac/pop_an
						background_af = background_ac/background_an
						popAFvsglobalAF.append((background_af, pop_af, background_ac, pop_ac, background_an, pop_an))
					vcf.read_line()
				
				print '%s LoF Allele Count: %s' % (args.pop, total_pop_ac)
				print '%s LoF Chromosome Count: %s' % (args.pop, total_pop_an)
				print 'Backround LoF Allele Count: %s' % total_background_ac
				print 'Background LoF Chromosome Count: %s' % total_background_an
				
				if args.output:
					with open(args.output, 'a') as out:
						out.write('BACKGROUND_AF\tPOP_AF\tBACKGROUND_AC\tPOP_AC\tBACKGROUND_AN\tPOP_AN\n')
						for line in popAFvsglobalAF:
							out.write('%s\n' % '\t'.join(map(str,line)))



def fishers_exact(acs, ans, idx):
	from scipy.stats import fisher_exact
	a = acs[idx]
	b = ans[idx]-a
	c = sum(acs)-a
	d = sum(ans)-c
	matrix = [[a, b], [c, d]]
	oddsratio, p_value = fisher_exact(matrix, alternative='greater')
	return p_value


def likelihood_ratio_test(acs, ans, idx, only_enriched=False):
	from math import log
	from scipy.stats import chi2
	def safe_log(x):
		if x == 0:
			return -1e9
		return log(x)
	global_af = float(sum(acs))/sum(ans)
	L_null = sum(acs)*safe_log(global_af)+(sum(ans)-sum(acs))*safe_log(1-global_af)
	ac_pop = acs[idx]
	an_pop = ans[idx]
	af_pop = float(ac_pop)/an_pop
	# time.sleep(0.5)
	if only_enriched and ac_pop < global_af*an_pop:
		return 1, 5
	ac_rest = sum(acs)-ac_pop
	an_rest = sum(ans)-an_pop
	af_rest = float(ac_rest)/an_rest
	# print af_pop, 1-af_pop, af_rest, 1-af_rest
	# get MLE log-likelihood
	L_alt = ac_pop*safe_log(af_pop)+(an_pop-ac_pop)*safe_log(1-af_pop)+ac_rest*safe_log(af_rest)+(an_rest-ac_rest)*safe_log(1-af_rest)
	chi_squared = 2*(L_alt-L_null)
	df = 1
	return 1 - chi2.cdf(chi_squared, df), chi_squared


def get_feature_from_annotation(annotation_list, key):
	this_set = get_set_from_annotation(annotation_list, key)
	if len(this_set) > 1:
		print >> sys.stderr, "WARNING: Set has more than one entry. Set was %s" % this_set
	return list(this_set)[0]


def filter_annotation_list(annotation_list, key, value_list, use_filter=True):
	if use_filter:
		return [x for x in annotation_list if x[key] in value_list]
	else:
		return [x for x in annotation_list if x[key] not in value_list]


def filter_annotation(annotation_list, key, value=None, use_filter=True):
	if value is None:
		if key.lower() == 'canonical':
			key = 'CANONICAL'
			value = 'YES'
		if key.lower() == 'lof':
			key = 'LoF'
			value = 'HC'
	if use_filter:
		return [x for x in annotation_list if x[key] == value]
	else:
		return [x for x in annotation_list if x[key] != value]


def get_set_from_annotation(annotation_list, key):
	return set([x[key] for x in annotation_list])


class Multinomial:
	def __init__(self, params):
		self.params = params
		self.combinations = []

	def exact_test(self, counts):
		n = sum(counts)
		k = len(counts)
		point = self.pmf(counts)
		self._permute(n, k, [])

		prob = point
		for case in self.combinations:
			p = self.pmf(case)
			if p < point:
				prob += p
		return prob

	def _permute(self, n, k, fixed):
		if k == 1:
			self.combinations.append(fixed + [n])
			return
		for a in range(0, n+1):
			self._permute(n-a, k-1, fixed + [a])


	def pmf(self, counts):
		if not(len(counts)==len(self.params)):
			raise ValueError("Dimensionality of count vector is incorrect")

		prob = 1.0
		for i,c in enumerate(counts):
			prob *= self.params[i]**counts[i]

		return prob * math.exp(self._log_multinomial_coeff(counts))

	def log_pmf(self,counts):
		if not(len(counts)==len(self._params)):
			raise ValueError("Dimensionality of count vector is incorrect")

		prob = 0.
		for i,c in enumerate(counts):
			prob += counts[i]*math.log(self._params[i])

		return prob + self._log_multinomial_coeff(counts)

	def _log_multinomial_coeff(self, counts):
		return self._log_factorial(sum(counts)) - sum(self._log_factorial(c)
														for c in counts)

	def _log_factorial(self, num):
		if not round(num)==num and num > 0:
			raise ValueError("Can only compute the factorial of positive ints")
		return sum(math.log(n) for n in range(1,num+1))



if __name__ == '__main__':
	import argparse
	import sys
	# sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools')
	# sys.path.insert(0, '/Users/birnbaum88/Desktop/Broad/ribosomal_occupancy/tools')
	from vcf_parsing import VCF_Parser
	from scipy.stats import chisquare
	import numpy as np
	from itertools import imap
	import math
	import time
	from collections import defaultdict

	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf', dest='vcf', help='VCF to parse.')
	parser.add_argument('--table', dest='table', help='Table containing Population AC/AN info.')
	parser.add_argument('--filter', dest='filter', action='store_true')
	parser.add_argument('--output', dest='output', required=False)
	parser.add_argument('--vcf_output', dest='vcf_output', required=False)
	parser.add_argument('--pop_spec', dest='pop_spec', action='store_true')
	parser.add_argument('--sig_thresh', dest='sig_thresh', type=float)
	parser.add_argument('--variant_level', dest='variant_level', action='store_true')
	parser.add_argument('--gene_level', dest='gene_level', action='store_true')
	parser.add_argument('--chisquare', dest='chisquare', action='store_true')
	parser.add_argument('--loo', dest='loo', action='store_true')
	parser.add_argument('--fisher', dest='fisher', action='store_true')
	parser.add_argument('--pop_max', dest='pop_max', action='store_true')
	parser.add_argument('--lrt', dest='lrt', action='store_true')
	parser.add_argument('--pop', dest='pop')
	parser.add_argument('--enriched', dest='enriched', action='store_true')
	args=parser.parse_args()
	main(args)






