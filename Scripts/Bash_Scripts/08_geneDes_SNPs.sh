#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=geneDes_SNPs
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 07/15/22                  Last Modified: 07/15/22 ###
###########################################################################
###########################################################################
###                     geneDes_SNPs.sh              					###
###########################################################################


cd $SLURM_SUBMIT_DIR
module load vcftools
module load bedops

# To generate distances between each SNP and its nearest genes upstream and downstream, you will need a .gff file 
# To prep the .gff file, first remove the trailing semicolons in the file and then remove all features that aren’t 
# genes (e.g., mRNA, exons, etc.)

cd /fs/scratch/PAS1533/smathur/wgr/geneDes/

awk '{gsub(/;$/,"");print}' Scate_HiC_rnd4.all.putative.function.gff > Scate_HiC.revisedgenome.all.gff
cat Scate_HiC.revisedgenome.all.gff | awk '$3 =="gene" {print $0}' > Scate_HiC.revisedgenome.all_justgenes.gff
cat Scate_HiC.revisedgenome.all.gff | awk '$3 =="CDS" {print $0}' > Scate_HiC.revisedgenome.all_justCDS.gff

# Use bedops to measure distances between genes and SNPs:

gff2bed < Scate_HiC.revisedgenome.all_justgenes.gff > Scate_HiC.revisedgenome.all_justgenes.bed
gff2bed < Scate_HiC.revisedgenome.all_justCDS.gff > Scate_HiC.revisedgenome.all_justCDS.bed
vcf2bed < ../vcf/onlyScat.10x.wholegenome.SNPs.recode.vcf > onlyScat.wgr.SNPs.bed
sort-bed onlyScat.wgr.SNPs.bed > onlyScat.wgr.sortedSNPs.bed
sort-bed Scate_HiC.revisedgenome.all_justgenes.bed > Scate_HiC.sortedgff_justgenes.bed
closest-features --delim '\t' --dist onlyScat.wgr.sortedSNPs.bed Scate_HiC.sortedgff_justgenes.bed > onlyScat.wgr.distance2genes_tab.txt


# Some markers will be upstream from a gene but there will be no genes downstream (or vice versa).  If we ignore this 
# we can end up incorporating a lot of SNPs that are at the extreme edge of a scaffold.  This makes me nervous because 
# for all we know there is a gene 100 bp downstream of the marker, but we can't tell because the scaffold doesn't extend 
# that far.  As such, I prefer to work with SNPs that have an annotated gene both upstream and downstream, so I know 
# exactly what kind of distances I'm working with.

# The following selects for SNPs that have a SNP both upstream and downstream 

cat onlyScat.wgr.distance2genes_tab.txt | awk '{if ($123!~/NA/) print $0;}' | awk '{if ($23!~/NA/) print $0;}' > distance2genes_noNA_clean.txt

# Keep only useful columns

cat distance2genes_noNA_clean.txt | cut -f1,2,3,122,123,124,125,131,132,133,134,135,136,142,143 > distance2genes_clean_all.txt

#######################################################################################################################################################################

# Rest of the downstream processing is done in R:

# Remove SNPs that are actually within annotated genes.
# Remove SNPs that are less than 1Mb on either side of an annotated gene


module load R
R
setwd("/fs/scratch/PAS1533/smathur/wgr/geneDes/") 
df <- read.csv("distance2genes_noGene_all.txt", header=F, sep="\t")
df.nogene <- df[-which(df$V9 == 0 | df$V15 == 0),] 

df.gendes <- df.nogene[which(abs(df.nogene$V9) > 1e6 & abs(df.nogene$V15) > 1e6), ] 
write.table(df.gendes, "Snps_in_genedeserts.txt",quote=F,row.names=F, col.names=F)

snps <- df.gendes[,c(1,2,3)]
write.table(snps, "genedeserts_snps.bed",quote=F,row.names=F, col.names=F)