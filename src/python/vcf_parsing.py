import re
import gzip
import argparse
import sys
import time
import os

class Thing:
    def __init__(self):
        self.x = 'yo'
    def __getitem__(self, y):
        return(self.x+'!'*y)
    
def try_convert_to_int(seq):
    ints = []
    for elem in seq:
        try:
            ints.append(int(elem))
        except ValueError:
            continue
    return ints

def generator(a_list):
    for elem in a_list:
        yield elem

class VCF_Parser:

    def __init__(self, vcf, subset=None):
        if type(vcf) == list:
            self.vcf = generator(vcf)
        else:
            self.vcf = gzip.open(vcf) if vcf.endswith('.gz') else open(vcf)
        self.output = open(subset, 'a') if subset else None
        vep_field_names = None
        header = None
        
        for line in self.vcf:
            self.line = line
            line = line.strip()

            # Reading header lines to get VEP and individual arrays
            if line.startswith('#'):
                if self.output:
                    self.output.write(self.line)
                line = line.lstrip('#')
                if line.find('ID=CSQ') > -1:
                    self.vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
                if line.startswith('CHROM'):
                    header = line.split('\t')
                    self.sample_names = header[9:]
                    self.header = dict(zip(header, range(len(header))))
                continue
            # Pull out annotation info from INFO and ALT fields
            self.fields = line.split('\t')
            self.info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', self.fields[self.header['INFO']])])
            self.ref_allele = self.fields[self.header['REF']]
            self.alt_allele = self.fields[self.header['ALT']]
            # Only reading lines with an annotation after this point
            if 'CSQ' not in self.info_field: 
                self.annotations = None
            else:
                self.annotations = [dict(zip(self.vep_field_names, x.split('|'))) for x in self.info_field['CSQ'].split(',')]
            break

        self.reading = True
    

    def __getitem__(self, col_name):
        return self.fields[self.header[col_name]]


    def read_line(self):
        try:
            self.line = self.vcf.next()
            self.fields = self.line.split('\t')
            self.ref_allele = self.fields[self.header['REF']]
            self.alt_allele = self.fields[self.header['ALT']]
            self.info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', self.fields[self.header['INFO']])])
            if 'CSQ' not in self.info_field: 
                self.annotations = None
            else:
                self.annotations = [dict(zip(self.vep_field_names, x.split('|'))) for x in self.info_field['CSQ'].split(',')]

        except StopIteration:
            self.reading = False


    def write_line(self):
        if self.output: 
            self.output.write(self.line)
            if not self.reading:
                self.output.close()


    def get_samples_carrying_variant(self, allele_num):
        mutants = []

        # Get index mapping for genotype fields  
        GT_format = self.fields[self.header['FORMAT']].split(':')
        GT_format = dict(zip(GT_format, range(len(GT_format))))

        # Identify mutants 
        for sample in self.sample_names:            
            gt_data = self.fields[self.header[sample]].split(':')
            gt = gt_data[GT_format['GT']]
            gt = try_convert_to_int(gt)
            if allele_num in gt:
                mutants.append((sample[2:], gt))
        return mutants


    def get_genotype_for_individual(self, indiv):
        GT_format = self.fields[self.header['FORMAT']].split(':')
        GT_format = dict(zip(GT_format, range(len(GT_format))))
        return self.fields[self.header[indiv]].split(':')[GT_format['GT']]


    def get_het_samples(self, gq_threshold=0):
        hets = []
        # Get index mapping for genotype fields  
        GT_format = self.fields[self.header['FORMAT']].split(':')
        GT_format = dict(zip(GT_format, range(len(GT_format))))

        for sample in self.sample_names:
            gt_data = self.fields[self.header[sample]].split(':')
            gt = gt_data[GT_format['GT']]
            gt = try_convert_to_int(gt)
            alt = filter(lambda x: x>0, gt)
            if 0 in gt and len(alt) > 0:
                gq = gt_data[GT_format['GQ']]
                if gq > gq_threshold:
                    hets.append((sample[2:], alt[0]-1))
        return hets


    def get_allele_by_num(self, allele_num):
        return self.fields[self.header['ALT']].split(',')[int(allele_num)-1]


    def get_variant(self):
        chrom = self.fields[self.header['CHROM']]
        pos = self.fields[self.header['POS']]
        return (chrom, int(pos))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', '--in', '-i', dest="vcf", help='Input VCF file (from VEP+LoF); may be gzipped', required=True)
    parser.add_argument('--var', dest="var", help='variant of interest', required=False)
    parser.add_argument('-o', '--rc', '--table', dest='output', help='Output file with gene rows and sample columns', required=False)
    parser.add_argument('--vcfsub', dest='vcf_sub', help='Output file of vcf lines corresponding to variants of interest', required=False)
    args = parser.parse_args()

    if args.var: 
        var = args.var.split(',')
        q_var = (var[0], int(var[1]))

    # with open('../genes/high_rpkm.genes', 'r') as f:
    #     expressed = [line.strip() for line in f]

    # bams = get_list_of_all_bams()
    # depths = get_depth_dict()

    vcf_parser = VCF_Parser(args.vcf)
    # vm = Variant_Mapper('protein_coding_genes.map')
    i = 0
    bed = []
    genes_being_queried = []
    vcf_lines = []

    while vcf_parser.reading:
        vcf_parser.read_line()

        # get variant
        var = vcf_parser.get_variant()
        if args.var:
            if q_var == var:
                # print vcf_parser.line
                samples = vcf_parser.get_samples_carrying_variant(1)
                samples = ['%s:%s' % (s[0], s[1]) for s in samples]
                print '\n'.join(samples)
                quit()
            else:
                continue

        # Only consider variants labeled as POSITIVE_TRAIN_SITE
        if 'POSITIVE_TRAIN_SITE' in vcf_parser.info_field:
            alleles = vcf_parser.get_field('ALT')
            # Out of these, only consider variants where there is ONE alternate allele
            if not len(alleles.split(',')) > 1:
                samples = vcf_parser.get_samples_carrying_variant(1)
                # Out of these, only consider variants where a SINGLE sample carries the mutation
                if len(samples) == 1:
                    genes = vm[var]
                    # Out of these, only consider variants that map to an annotated gene
                    if genes:
                        overlap = filter(lambda gene: gene['id'] in expressed, genes)
                        # Out of these, only consider variants that map to annotated genes with high average RPKM
                        if overlap:
                            # At this point, we have filtered down to 327 variants
                            regions = [(g['chrom'], g['start'], g['end']) for g in overlap]
                            gene_names = [g['id'] for g in overlap]
                            bed += regions
                            genes_being_queried += gene_names
                            vcf_lines.append(vcf_parser.get_line())


    # print "Querying %s regions." % len(bed) 
    # results = parallelize_coverageBed(bed, bams)
    # samples = [bam.split('/')[-1].split('_')[0] for bam in bams]

    if args.output:
        from profiling import write_read_count_table
        # write table to file
        write_read_count_table(args.output, results, bams, genes_being_queried)

    if args.vcf_sub:
        with open(args.vcf_sub, 'a') as info:
            for line in vcf_lines:
                info.write(line)








