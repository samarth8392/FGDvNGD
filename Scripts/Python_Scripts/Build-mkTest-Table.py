#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from fuc import pyvcf
import pandas as pd
from copy import deepcopy 

parser = argparse.ArgumentParser(description='Builds McDonald-Kreitman contingency tables for detecting selection in coding sequences')
parser.add_argument("-s","--seqfile",
					type=str,
					default='',
					help="Sequence file in fasta format")
parser.add_argument("-g","--gff",
					type=str,
					default='',
					help="gff3 file containing genes as a type")
parser.add_argument("-v","--vcf",
					type=str,
					help="VCF file with variant calls")
parser.add_argument("-ig","--ingroup",
					type=str,
					help="text file with ingroup samples")
parser.add_argument("-og","--outgroup",
					type=str,
					help="text file with outgroup samples")
parser.add_argument("-o","--output",
					type=str,
					default='',
					help="name of output csv file")
args=parser.parse_args()

seqfile = args.seqfile
gff = args.gff
vcf = args.vcf
feature = "CDS"
output = args.output
ingroup = args.ingroup
outgroup = args.outgroup

## This is a dictionary of Tsil and Trepl scores for SnIPRE. Tsil will be first, Trepl second
AA_scores = {'A':[1,2], 'R':[1,2], 'N':[0.333,2.667], 'D':[0.333,2.667], 'C':[0.333,2.667], 'E':[0.333,2.667], 'Q':[0.333,2.667], 'G':[1,2], 'H':[0.333,2.667], 'I':[0.667,2.333], 'L':[1.333,1.667], 'K':[0.333,2.667], 'M':[0,3], 'F':[0.333,2.667], 'P':[1,2], 'S':[1,2], 'T':[1,2], 'W':[0,3], 'Y':[0.333,2.667], 'V':[1,2], '*':[0.667,2.333]}

#########################################################################################

class GFF:
    def __init__(self, gff_entry):
        self.seqid = gff_entry[0]
        self.source = gff_entry[1]
        self.type = gff_entry[2]
        self.start = int(gff_entry[3])
        self.end = int(gff_entry[4])
        self.score = gff_entry[5]
        self.strand = gff_entry[6]
        self.phase = gff_entry[7]
        if ';' in gff_entry[8]:
            self.attributes = parse_GFF_attributes(gff_entry[8])
        else:
            self.attributes = gff_entry[8]


class gff_attributes:
    def __init__(self, gff_attributes):
        self.ID = gff_attributes[0]
        self.Name = gff_attributes[1]
        self.Alias = gff_attributes[2]
        self.Note = gff_attributes[3]



class CDS_object:
    def __init__(self,CDS_object):
        self.ID = CDS_object[0]
        self.seq = CDS_object[1]
        self.exons = CDS_object[2]
        self.strand = CDS_object[3]


def GFF_parse(gff) :
	gff_list = []
	with open(gff) as OF:
		reader = csv.reader(OF, delimiter='\t')
		for row in reader :
			if (row[0][0] != '#') and (len(row) > 1) :
				gff_entry = GFF(row)
				gff_list.append(gff_entry)
	return(gff_list)


def parse_GFF_attributes(gff_attribute) :
    split_attributes_1 = list(gff_attribute.split(';'))
    split_attributes = []
    for attribute in split_attributes_1:
        split_attributes.append(list(attribute.split('=')))
    ID = [x[1] for x in split_attributes if x[0] == 'ID']
    Name = [x[1] for x in split_attributes if x[0] == 'Name']
    Alias = [x[1] for x in split_attributes if x[0] == 'Alias']
    Note = [x[1] for x in split_attributes if x[0] == 'Note']
    attributes_list = gff_attributes([ID, Name, Alias, Note])
    return(attributes_list)


def Pull_CDS(sequence_list, gff_list) :
    gene_gff_list = [ entry for entry in gff_list if entry.type == 'gene']
    CDS_gff_list = [ entry for entry in gff_list if entry.type == 'CDS']       
    CDS_obj_dic = {}
    CDS_ID = ''
    for entry in CDS_gff_list:
        if entry.attributes.ID[0] == CDS_ID:
            exon_start = entry.start -1
            exon_end = entry.end
            new_exon = SeqRecord(
                seq[0].seq[exon_start:exon_end],
                id=entry.attributes.ID[0])
            if entry.strand == '-':
                rc_seq = new_exon.seq.reverse_complement()
                new_exon.seq = rc_seq
            extended_CDS = new_seq.seq + new_exon.seq
            new_seq.seq = extended_CDS
            exon = [entry.start,entry.end]
            exons.append(exon)
            CDS_entry = CDS_object([CDS_id,new_seq,exons,entry.strand])
        else:
            if entry != CDS_gff_list[0]:
                CDS_obj_dic[CDS_id] = CDS_entry
            CDS_id = entry.attributes.ID[0].split('-RA')[0]
            seq = [seq for seq in sequence_list if seq.id == entry.seqid]
            gene = [gene for gene in gene_gff_list if CDS_id == gene.attributes.ID[0]]
            if len(seq) > 1:
                print('so, somehow there is a gff that matches more than one sequence')
                print(entry.attributes.ID)
                break
            CDS_ID = entry.attributes.ID[0]
            exon_start = entry.start -1
            exon_end = entry.end
            if len(gene[0].attributes.Note) > 0 :
                new_seq = SeqRecord(
                    seq[0].seq[exon_start:exon_end],
                    id=entry.attributes.ID[0],
                    description=gene[0].attributes.Note[0]
                )
            else:
                new_seq = SeqRecord(
                    seq[0].seq[exon_start:exon_end],
                    id=entry.attributes.ID[0]
                )
            if entry.strand == '-':
                rc_seq = new_seq.seq.reverse_complement()
                new_seq.seq = rc_seq
            exons = [[entry.start,entry.end]]
            CDS_entry = CDS_object([CDS_id,new_seq,exons,entry.strand])
    CDS_obj_dic[CDS_id] = CDS_entry
    return(CDS_obj_dic)


def genotype_from_gt_string(gt_string):
    if '.' in gt_string:
        haplotype_1 = '.'
        haplotype_2 = '.'
    elif '/' in gt_string:
        haplotype_1 = gt_string.split('/')[0] 
        haplotype_2 = gt_string.split('/')[1]
    elif '|' in gt_string:
        haplotype_1 = gt_string.split('|')[0] 
        haplotype_2 = gt_string.split('|')[1]
    else:
        print(gt_string)  
    if haplotype_1 != '.':
            haplotype_1 = int(haplotype_1)
    if haplotype_2 != '.':
            haplotype_2 = int(haplotype_2)          
    gt = [haplotype_1, haplotype_2]
    return(gt)



def fixed_in_to_out(in_gt,out_gt):
    if ['.', '.'] in in_gt:
        in_gt = [ x for x in in_gt if x != ['.', '.'] ]
    uniq_in_gt = [list(x) for x in set(tuple(x) for x in in_gt)]
    if (uniq_in_gt == [[1,1]]) or (uniq_in_gt == [[0,0]]):
        return(True)
    else:
        return(False)



def CDS_pos_from_ref_position(snp_position, CDS_object):
    CDS_seq_pos = 0
    exons_length = [exon[1]+1-exon[0] for exon in CDS_object.exons]
    exon_counter = 0
    exon_addition = 0
    for exon in CDS_object.exons:
        if exon[0] <= snp_position <= exon[1]:
            if CDS_object.strand == '+':
                exon_addition = snp_position - exon[0]
            elif CDS_object.strand == '-':
                 exon_addition = (exon[1] - snp_position)
            else :
                print('uhhh, this doesnt have strand info?')
            break
        elif exon == CDS_object.exons[-1] and (not exon[0] <= snp_position <= exon[1]):
            print('Hmmm. Looks like this snp isnt in an exon')
            break
        else:
            exon_counter += 1
    for i in range(0,exon_counter):
        CDS_seq_pos += exons_length[i]
    CDS_seq_pos += exon_addition
    return(CDS_seq_pos)


def test_synonymous_nonsynonymous(snp_vcf_entry, CDS_snp_pos, CDS_object):
    reference_AA_sequence = deepcopy(CDS_object.seq.seq.translate())
    test_seq = deepcopy(CDS_object.seq)
    if CDS_object.strand == '+':
        replacement_seq = test_seq.seq[:CDS_snp_pos] + snp_vcf_entry['ALT'] + test_seq.seq[CDS_snp_pos+1:]
    elif CDS_object.strand == '-':
        if snp_vcf_entry['ALT'] == 'A':
            replacement_seq = test_seq.seq[:CDS_snp_pos] + 'T' + test_seq.seq[CDS_snp_pos+1:]
        elif snp_vcf_entry['ALT'] == 'T':
            replacement_seq = test_seq.seq[:CDS_snp_pos] + 'A' + test_seq.seq[CDS_snp_pos+1:]
        elif snp_vcf_entry['ALT'] == 'C':
            replacement_seq = test_seq.seq[:CDS_snp_pos] + 'G' + test_seq.seq[CDS_snp_pos+1:]
        elif snp_vcf_entry['ALT'] == 'G':
            replacement_seq = test_seq.seq[:CDS_snp_pos] + 'C' + test_seq.seq[CDS_snp_pos+1:]
    test_seq.seq = replacement_seq
    test_AA_sequence = test_seq.seq.translate()
    if reference_AA_sequence == test_AA_sequence:
        substitution = 'synonymous'
    elif reference_AA_sequence != test_AA_sequence:
        substitution = 'nonsynonymous'
    return(substitution)


def calc_Tsil(CDS_object):
    AA_seq = CDS_object.seq.seq.translate()
    Tsil = 0
    for AA in AA_seq:
        Tsil += AA_scores[AA][0]
    return(Tsil)


def calc_Trepl(CDS_object):
    AA_seq = CDS_object.seq.seq.translate()
    Trepl = 0
    for AA in AA_seq:
        Trepl += AA_scores[AA][1]
    return(Trepl)

########################################################################################

print('Reading sample lists')
with open(ingroup) as ingroup_file:
    ingroup_samples = ingroup_file.readlines()

ingroup_samples = [ x.strip() for x in ingroup_samples ]

len_ingroup = 2*(len(ingroup_samples))

with open(outgroup) as outgroup_file:
    outgroup_samples = outgroup_file.readlines()

outgroup_samples = [ x.strip() for x in outgroup_samples ]

len_outgroup = 2*(len(outgroup_samples))

# Combine ingroup and outgroup for filtering
all_samples = ingroup_samples + outgroup_samples
print(f'Ingroup samples: {len(ingroup_samples)}')
print(f'Outgroup samples: {len(outgroup_samples)}')
print(f'Total samples to keep: {len(all_samples)}')

print('reading gff file')
## 1) Read in GFF. Reduce to CDS annotations
gff_list = GFF_parse(gff)
CDS_gff_list = [ entry for entry in gff_list if entry.type == 'CDS']

print('reading vcf file')
## 2) Read in VCF
vcf_df = pyvcf.VcfFrame.from_file(vcf)

#######################################################################################
# OPTIMIZATION 0: Filter VCF to only include ingroup and outgroup samples
# This reduces memory usage and speeds up all downstream operations
#######################################################################################
print('Filtering VCF to only ingroup and outgroup samples (OPTIMIZATION)')
print(f'Original VCF columns: {len(vcf_df.df.columns)}')

# Get list of all sample columns in VCF (excluding standard VCF columns)
vcf_standard_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
all_vcf_samples = [col for col in vcf_df.df.columns if col not in vcf_standard_cols]

# Find which samples to keep (intersection of VCF samples and our ingroup/outgroup)
samples_to_keep = [s for s in all_samples if s in all_vcf_samples]
missing_samples = [s for s in all_samples if s not in all_vcf_samples]

if missing_samples:
    print(f'WARNING: {len(missing_samples)} samples not found in VCF: {missing_samples[:5]}...')
print(f'Keeping {len(samples_to_keep)} samples from VCF')

# Filter VCF to only these columns
cols_to_keep = vcf_standard_cols + samples_to_keep
vcf_df.df = vcf_df.df[cols_to_keep]
print(f'Filtered VCF columns: {len(vcf_df.df.columns)}')

print('reading sequences')
## I think having a chromosome dictionary might be useful
sequences = list(SeqIO.parse(seqfile,"fasta"))
seq_dict = {}
for seq in sequences:
    seq_dict[seq.id] = seq.seq

print('finding unique CDSs')
uniq_CDSs = []
for CDS in CDS_gff_list:
    if CDS.attributes.ID[0].split('-RA')[0] not in uniq_CDSs:
        uniq_CDSs.append(CDS.attributes.ID[0].split('-RA')[0])


print('making CDS dictionary')
CDS_dic = {}
for CDS in uniq_CDSs:
    CDS_parts = [ x for x in CDS_gff_list if x.attributes.ID[0].split('-RA')[0] == CDS]
    CDS_dic[CDS] = CDS_parts

#######################################################################################
# OPTIMIZATION 1: Pre-filter VCF by chromosome and cache
# This avoids deep copying the entire VCF for every exon
#######################################################################################
print('pre-filtering VCF by chromosome (OPTIMIZATION)')
chrom_vcf_cache = {}
unique_chroms = set([exon.seqid for CDS in uniq_CDSs for exon in CDS_dic[CDS]])
for chrom in unique_chroms:
    chrom_vcf_cache[chrom] = vcf_df.df[vcf_df.df['CHROM'] == chrom].copy()
    print(f'  Cached {len(chrom_vcf_cache[chrom])} SNPs for {chrom}')

#######################################################################################
# OPTIMIZATION 2: Use cached chromosome VCFs instead of deep copying entire VCF
#######################################################################################
CDS_snp_counter = 0
print('making CDS snp dataframes (OPTIMIZED)')
CDS_snp_dataframes = {}
for CDS in uniq_CDSs:
    CDS_snp_counter += 1
    if CDS_snp_counter % 100 == 0:  # Progress update every 100 CDSs
        print('on CDS ', CDS, ' which is ', CDS_snp_counter, ' of ', len(uniq_CDSs))
    snps_list = []
    for exon in CDS_dic[CDS]:
        # USE CACHED CHROMOSOME VCF instead of deepcopy(vcf_df.df)!
        chrom_snps = chrom_vcf_cache[exon.seqid]
        exon_snps = chrom_snps[(chrom_snps['POS'] >= exon.start) & (chrom_snps['POS'] <= exon.end)]
        snps_list.append(exon_snps)
    CDS_snp_dataframes[CDS] = snps_list

print('making CDS objects')
CDS_seqs = Pull_CDS(sequences, gff_list)

SnIPRE_table = [["geneID", "PR", "FR", "PS", "FS", "Tsil", "Trepl", "nout","npop"]]

#######################################################################################
# OPTIMIZATION 3: Use itertuples() instead of iterrows() for faster iteration
#######################################################################################
print('starting calculations')
gene_counter = 0
for CDS in uniq_CDSs:
    gene_counter += 1
    if gene_counter % 100 == 0:  # Progress update
        print(f'Processing gene {gene_counter} of {len(uniq_CDSs)}: {CDS}')
    
    nonsynonymous_fixed = 0
    nonsynonymous_polymorphic = 0
    synonymous_fixed = 0
    synonymous_polymorphic = 0
    
    for exon_snps in CDS_snp_dataframes[CDS]:
        # OPTIMIZATION: Use itertuples() instead of iterrows()
        for snp in exon_snps.itertuples(index=False):
            snp_dict = snp._asdict()  # Convert to dictionary for easier access
            
            out_gt = []
            for sample in outgroup_samples:
                if sample in snp_dict:  # Check if sample exists after filtering
                    gt_string = snp_dict[sample].split(':')[0]
                    gt = genotype_from_gt_string(gt_string)
                    out_gt.append(gt)
            
            in_gt = []
            for sample in ingroup_samples:
                if sample in snp_dict:  # Check if sample exists after filtering
                    gt_string = snp_dict[sample].split(':')[0]
                    gt = genotype_from_gt_string(gt_string)
                    in_gt.append(gt)
            
            if fixed_in_to_out(in_gt, out_gt):
                fixed_polymorphic = 'fixed'
            else:
                fixed_polymorphic = 'polymorphic'
            
            snp_pos = CDS_pos_from_ref_position(snp_dict['POS'], CDS_seqs[CDS])
            synon_nonsynon = test_synonymous_nonsynonymous(snp_dict, snp_pos, CDS_seqs[CDS])
            
            if fixed_polymorphic == 'fixed' and synon_nonsynon == 'synonymous':
                synonymous_fixed += 1
            elif fixed_polymorphic == 'fixed' and synon_nonsynon == 'nonsynonymous':
                nonsynonymous_fixed += 1
            elif fixed_polymorphic == 'polymorphic' and synon_nonsynon == 'synonymous':
                synonymous_polymorphic += 1
            elif fixed_polymorphic == 'polymorphic' and synon_nonsynon == 'nonsynonymous':
                nonsynonymous_polymorphic += 1
    
    Tsil = calc_Tsil(CDS_seqs[CDS])
    Trepl = calc_Trepl(CDS_seqs[CDS])
    entry = [CDS_seqs[CDS].ID, nonsynonymous_polymorphic, nonsynonymous_fixed, 
             synonymous_polymorphic, synonymous_fixed, Tsil, Trepl, len_outgroup, len_ingroup]
    print(entry)
    SnIPRE_table.append(entry)

with open(output, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    for row in SnIPRE_table:
        csv_writer.writerow(row)

csv_file.close()

print(f'\nDone! Output written to {output}')