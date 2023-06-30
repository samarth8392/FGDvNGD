#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH --job-name=ldNe
#SBATCH -e ldNe
#SBATCH -o ldNe

##########################################################################
###                          Samarth Mathur, PhD                        ###
###                       The Ohio State University                     ###
###                                                                     ###
###     Date Created: 03/30/22                  Last Modified: 08/22/22 ###
###########################################################################
###########################################################################
###                     ldNe.sh                        					###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load vcftools
#module load plink/1.90b6.4

### Comment:
### Adding Waples contemporary Ne estimates would be very useful and help support the argument (i.e. line 524-525). 
### His 2016 paper adds an adjustment factor for genome data BUT I've written a script that simply randomly samples 
### the genome X number of times to calculate Ne (via R2) with SE. It works well and is quick and I'm happy to share via the editor, 
### as having pi (or long-term Ne) and contemporary Ne I think would help round-off the ROH patterns. That said, it's not difficult to write: 
### Waples (2006) equations are based on genotype correlations (R2) and you can apply the sample size correction of Waples et al. (2016). 
### All we did was randomly sample the VCF file to produce "unlinked" data of a desired size.



# Step 1: Randomly choose 13000 SNPs from gene desert SNPs (see script 08_geneDes_SNPs) and calculate LD for each population and repeat it 100 times


for i in {1..100}
do
	shuf -n 13000 genedeserts_snps.bed > randsites/genDes.rand.run${i}.sites

	for pop in BPNP CCRO CEBO JENN KBPP KLDR MOSQ PRDF ROME SPVY SSSP WLRD
	do
		vcftools --gzvcf all121.final.10x.final.SNPs.vcf.gz \
		--maf 0.05 --geno-r2 \
		--keep ${pop}.all121.sample.list \
		--positions randsites/genDes.rand.run${i}.sites \
		--out ${pop}/all121.${pop}.run${i}
	done
done


# Step 2: Get LDNe using Rscript LDNe.R
