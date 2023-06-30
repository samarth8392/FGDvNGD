#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH --job-name=snpeff_all121
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 08/18/21                  Last Modified: 09/10/21 ###
###########################################################################
###########################################################################
###                     variant_annotation.sh                        	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load snpeff/4.2

# SNPeff

# Steps:
# 1. Copy snpeff.config file from /group/bioinfo/apps/apps/snpEff-4.3/ to local directory
# 2. Make a directory in your local directory called "data" and another directory within data with your species name. I call it scat (S.catanatus)
# 3. copy ref fasta and annotation gff as sequence.fa and genes.gff into "scat" folder
# 4. Change the following in config file:
# 	(a) Add in the first two lines:
#			# S.catanatus genome
#			scat.genome : S.catanatus
#	(b) Change data.dir to /fs/scratch/PAS1533/smathur/wgr/variant_annotation/snpeff/data/
# 5. Build the database

cd /fs/scratch/PAS1533/smathur/wgr/variant_annotation/snpeff/
java -jar /usr/local/snpeff/4.2/snpEff/snpEff.jar build \
-c snpEff.config -gff3 -v scat &> build.logfile.txt

# If the database builds without error, you should see snpEffectPredictor.bin inside your scat folder

# Step6: Annotate your variants

java -jar /usr/local/snpeff/4.2/snpEff/snpEff.jar ann \
-stats -c snpEff.config \
-no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v scat \
/fs/scratch/PAS1533/smathur/wgr/vcf/all121.10x.wholegenome.SNPs.recode.vcf \
> SNPeff.all121.final.10x.vcf


## Filter different sites

cat SNPeff.all121.final.10x.vcf | \
java -jar /usr/local/snpeff/4.2/snpEff/SnpSift.jar filter "( EFF[*].EFFECT = 'missense_variant' )" | \
grep -v "WARNING" > SNPeff.all121.final.nonSyn.vcf

cat SNPeff.all121.final.10x.vcf | \
java -jar /usr/local/snpeff/4.2/snpEff/SnpSift.jar filter "( EFF[*].IMPACT = 'HIGH' )" | \
grep -v "WARNING" > SNPeff.all121.final.LoF.vcf

cat SNPeff.all121.final.10x.vcf | \
java -jar /usr/local/snpeff/4.2/snpEff/SnpSift.jar filter "( EFF[*].EFFECT = 'synonymous_variant' )" | \
grep -v "WARNING" > SNPeff.all121.final.Syn.vcf
