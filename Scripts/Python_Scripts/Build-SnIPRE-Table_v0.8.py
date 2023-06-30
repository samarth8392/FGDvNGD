#!/usr/bin/env python


import argparse
import csv
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from fuc import pyvcf
import pandas as pd
from copy import deepcopy 

parser = argparse.ArgumentParser(description='writes a table of CG content in a genome (although technically it will work with any fasta)')
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
					help="feature you would like to collect in a fasta file. Must be a string matched in your gff")
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
    #feature_seqs = []
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
                #feature_seqs.append(new_seq)
                #feature_seqs.append(CDS_object)
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
    #feature_seqs.append(new_seq)
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
        #print('missing data in ingroup. will be removed')
        in_gt = [ x for x in in_gt if x != ['.', '.'] ]
    #if ['.', '.'] in out_gt:
        #print('missing data in ingroup. will be removed')
    #    out_gt = [ x for x in out_gt if x != ['.', '.'] ]
    uniq_in_gt = [list(x) for x in set(tuple(x) for x in in_gt)]
    #print(uniq_in_gt)
    #uniq_out_gt = [list(x) for x in set(tuple(x) for x in out_gt)]
    ## Note, here we can make the assumption that if an allele is fixed in the ingroup, it is fix and derived even if the ingroup allele is the reference
    ## The reason being we are assuming an input dataset of biallelic snps, so if a snp is represented in the
    ## dataframe then it must be variable, and if it isn't variable in the ingroup, it must be variable in the outgroup
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
    #print(CDS_object.seq.seq[CDS_snp_pos], CDS_object.strand)
    #print(snp_vcf_entry['ALT'])
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
    #print(CDS_object.seq.seq == replacement_seq)
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

with open(ingroup) as ingroup_file:
    ingroup_samples = ingroup_file.readlines()

ingroup_samples = [ x.split('\n')[0] for x in ingroup_samples ]

len_ingroup = 2*(len(ingroup_samples))

with open(outgroup) as outgroup_file:
    outgroup_samples = outgroup_file.readlines()

outgroup_samples = [ x.split('\n')[0] for x in outgroup_samples ]

len_outgroup = 2*(len(outgroup_samples))

print('reading gff file')
## 1) Read in GFF. Reduce to CDS annotations
gff_list = GFF_parse(gff)
CDS_gff_list = [ entry for entry in gff_list if entry.type == 'CDS']

print('reading vcf file')
## 2) Read in VCF. Reduce to SNPs that occur in CDS regions
vcf_df = pyvcf.VcfFrame.from_file(vcf)

#print('reducing vcf to coding snps')
#coding_snps = pd.DataFrame(columns = vcf_df.df.columns)
#seq_name = ''
#for x in range(0,len(vcf_df.df)):
#for x in range(0,10):
#    print(x, ' of ', str(len(vcf_df.df)), ' snps completed')
    #print(vcf_df.df.loc[x]['CHROM'])
#    if vcf_df.df.loc[x]['CHROM'] != seq_name:
#        seq_name = vcf_df.df.loc[x]['CHROM']
#        CDSs = [ CDS for CDS in CDS_gff_list if CDS.seqid == seq_name]
        #print(len(CDSs))
#    CDSs = [ CDS for CDS in CDSs if CDS.end >= vcf_df.df.loc[x]['POS']]
   #print(len(CDSs))
#    for CDS in CDSs:
#        if vcf_df.df.loc[x]['POS'] >= CDS.start:
            #vcf_df.df.loc[x]
#            coding_snps = pd.concat([coding_snps,vcf_df.df.loc[x].to_frame().T], ignore_index = True)
            #coding_snps.append(vcf_df.df.loc[x].to_frame().T)
#            break

print('reading sequences')
## I think having a chromosome dictionary might be useful
sequences = list(SeqIO.parse(seqfile,"fasta"))
seq_dict = {}
for seq in sequences:
    seq_dict[seq.id] = seq.seq


## 3) For each gene, make a list of variants, determine their codon, determine if alternate allele is synonmous or nonsynonomous
## 3a) also determine if alleles fixed or polymoprhic.
## 3b) calculate the Tsil and Trepl values


## I think here are the steps
## 1) make a table or dictionary or whatever of CDS

## 2) Convert the SNP information into a 'CDS' based vcf where positions denote positions in the CDS.

## 3) Identify the codon for each SNP and whether it is synonymous or nonsynonymous

## 4) Identify whether the snp is fixed or not
print('finding unique CDSs')
uniq_CDSs = []
for CDS in CDS_gff_list:
    #print(CDS.attributes.ID[0].split('-RA')[0])
    if CDS.attributes.ID[0].split('-RA')[0] not in uniq_CDSs:
        uniq_CDSs.append(CDS.attributes.ID[0].split('-RA')[0])


print('making CDS dictionary')
CDS_dic = {}
for CDS in uniq_CDSs:
    CDS_parts = [ x for x in CDS_gff_list if x.attributes.ID[0].split('-RA')[0] == CDS]
    CDS_dic[CDS] = CDS_parts


#CDS_snp_dataframes = {}
#for CDS in uniq_CDSs:
#    snps_list = []
#    for exon in CDS_dic[CDS]:
#        exon_snps = coding_snps[coding_snps['POS'] >= exon.start]
#        exon_snps = exon_snps[ exon_snps['POS'] <= exon.end ]
#        snps_list.append(exon_snps)
#    CDS_snp_dataframes[CDS] = snps_list

CDS_snp_counter = 0
print('making CDS snp dataframes')
CDS_snp_dataframes = {}
for CDS in uniq_CDSs:
    CDS_snp_counter += 1
    print('on CDS ', CDS, ' which is ', CDS_snp_counter, ' of ', len(uniq_CDSs))
    snps_list = []
    for exon in CDS_dic[CDS]:
        tmp_snps_list = deepcopy(vcf_df.df)
        exon_snps = tmp_snps_list[tmp_snps_list['CHROM'] == exon.seqid]
        exon_snps = exon_snps[exon_snps['POS'] >= exon.start]
        exon_snps = exon_snps[ exon_snps['POS'] <= exon.end ]
        snps_list.append(exon_snps)
    CDS_snp_dataframes[CDS] = snps_list

print('making CDS objects')
CDS_seqs = Pull_CDS(sequences, gff_list)

## So, now I have two dictionaries with 'CDS' names as keys.
## One links all the sequence and exon information to the CDS, the other links the SNP info to the CDS.

SnIPRE_table = [["geneID", "PR", "FR", "PS", "FS", "Tsil", "Trepl", "nout","npop"]]

print('starting calculations')
for CDS in uniq_CDSs:
#for CDS in test_uniq_CDSs:
    nonsynonymous_fixed = 0
    nonsynonymous_polymorphic = 0
    synonymous_fixed = 0
    synonymous_polymorphic = 0
    for exon_snps in CDS_snp_dataframes[CDS]:
        for snp in exon_snps.iterrows():
            #print(snp[1]['POS'])
            out_gt = []
            for sample in outgroup_samples:
                gt_string = snp[1][sample].split(':')[0]
                gt = genotype_from_gt_string(gt_string)
                out_gt.append(gt)
            in_gt = []
            for sample in ingroup_samples:
                gt_string = snp[1][sample].split(':')[0]
                gt = genotype_from_gt_string(gt_string)
                in_gt.append(gt)
            #print(set(out_gt))
            #print(set(in_gt))
            if fixed_in_to_out(in_gt, out_gt) :
                fixed_polymorphic = 'fixed'
            else :
                fixed_polymorphic = 'polymorphic'
            #print(fixed_polymorphic)
            snp_pos = CDS_pos_from_ref_position(snp[1]['POS'], CDS_seqs[CDS])
            synon_nonsynon = test_synonymous_nonsynonymous(snp[1], snp_pos, CDS_seqs[CDS])
            #print(synon_nonsynon)
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
    entry = [CDS_seqs[CDS].ID, nonsynonymous_polymorphic, nonsynonymous_fixed, synonymous_polymorphic, synonymous_fixed, Tsil, Trepl, len_outgroup, len_ingroup]
    #print(Tsil+Trepl)
    #print(len(CDS_seqs[CDS].seq.seq))
    print(entry)
    SnIPRE_table.append(entry)
    
        

with open(output, 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter = ',')
    for row in SnIPRE_table:
        csv_writer.writerow(row)


csv_file.close()
