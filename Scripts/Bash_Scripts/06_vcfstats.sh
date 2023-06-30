#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH --job-name=vcfstats
#SBATCH -e vcfstats
#SBATCH -o vcfstats

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/02/21                  Last Modified: 02/21/22 ###
###########################################################################
###########################################################################
###                     vcfstats.sh              						###
###########################################################################


cd $SLURM_SUBMIT_DIR
module load vcftools

cd /fs/ess/scratch/PAS1533/smathur/wgr/vcfstats/

##### Variant statistics ####

### Whole genome 

for i in depth site-pi het relatedness2 missing-indv missing-site
do
	vcftools --gzvcf all121.final.10x.final.SNPs.vcf.gz \
	--$i \
	--out stats/all121.allchr
done

# Get vcfs by pop

cd /fs/ess/scratch/PAS1533/smathur/wgr/vcf/

while read -a pop
do
	vcftools --gzvcf all121.final.10x.final.SNPs.vcf.gz \
	--recode --recode-INFO-all \
	--keep /fs/ess/scratch/PAS1533/smathur/wgr/lists/bypop/${pop}.all121.sample.list \
	--out bypop/${pop}.all121.allchr.finalSNPs
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames.txt