#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 150:00:00
#SBATCH --job-name=dos
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 11/08/21                  Last Modified: 12/16/21 ###
###########################################################################
###########################################################################
###                     DOS.sh                        				    ###
###########################################################################

cd $SLURM_SUBMIT_DIR

MAINDIR="/fs/ess/scratch/PAS1533/smathur/Scat_Popgen"
cd $MAINDIR

# module load vcftools

# Calculating DOS values using variants only in KLDR population (and STER as outgroup) & removing sites that are monomorphic ˛

# Step1: Get vcf from KLDR and STER and remove monomorphs

vcftools --gzvcf all121.final.10x.final.SNPs.vcf.gz \
--keep scat_pop_KLDR.list \
--keep ster_pop_STER.list \
--recode --recode-INFO-all \
--out onlyKLDR

vcftools --vcf onlyKLDR.recode.vcf \
--freq \
--out onlyKLDR


vcftools --vcf onlyKLDR.recode.vcf \
--recode --recode-INFO-all \
--exclude-positions monomorphs.txt \
--out onlyKLDR_nofixed


# Step2: Split vcf and gff by sacffold

INPUT_GFF="Scate_HiC_rnd4.all.putative.function.gff"

cd $MAINDIR/scaffoldGFFs
cut -f1 "$INPUT_GFF" | sort -u | while read scaffold; do
    grep "^${scaffold}" "$INPUT_GFF" > "${scaffold}.gff"
done

cd $MAINDIR/scaffoldGFFs/
for chr in $(ls *.gff | sed -r 's/.gff//'  | uniq)
do
	cd $MAINDIR/ScaffoldVCFs

	vcftools --vcf onlyKLDR_nofixed.recode.vcf \
	--chr $chr --recode --recode-INFO-all \
	--out onlyKLDR.$chr
done

# Step3: Create and activate conda environment

# Create environment from YAML file
conda env create -f environment.yml

# Activate the environment
source ~/.bashrc
conda activate mkTest_env

#Step3: Run python code written by Andrew Mason to get continegency table for DOS calculation.

# -s Sequence file in fasta format
# -g gff3 file containing genes as a type
# -v feature you would like to collect in a fasta file. Must be a string matched in your gff
# -ig text file with ingroup samples
# -og text file with outgroup samples
# -o name of csv output file

cd $MAINDIR/scaffoldGFFs/
for chr in $(ls *.gff | sed -r 's/.gff//'  | uniq)
do
    cd $MAINDIR/tables

    python Build-mkTest-Table.py \
    -s ../Scatenatus_HiC_v1.1.fasta \
    -g $MAINDIR/scaffoldGFFs/${chr}.gff \
    -v $MAINDIR/ScaffoldVCFs/${chr}.vcf \
    -ig ../kldr.txt \
    -og ../outgroup.txt \
    -o onlyKLDR_${chr}_2x2.csv
done

cd $MAINDIR/tables
cat *.csv > KLDR.ig_STER.og.2x2table.txt


# Next run R_scripts/DoS.R to calculate DoS values and plot the results.