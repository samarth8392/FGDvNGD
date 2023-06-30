#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=base_recal
#SBATCH -e %x_%j
#SBATCH -o %x_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 02/02/21                  Last Modified: 07/27/21 ###
###########################################################################
###########################################################################
###                     base_recal.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module load bwa
module load samtools/1.8
module load picard/2.18.17
module load gatk #gatk/4.1.2.0
module load vcftools 

# Base recalibration
#The base quality score recalibrator treats every reference mismatch as indicative of machine error. 
#True polymorphisms are legitimate mismatches to the reference and shouldn't be counted against the quality 
#of a base. We use a database of known polymorphisms to skip over most polymorphic sites. 
#Unfortunately without this information the data becomes almost completely unusable since the quality of 
#the bases will be inferred to be much much lower than it actually is as a result of the reference-mismatching SNP sites. 

#However, you can bootstrap a database of known SNPs. Here's how it works: 
#First do an initial round of SNP calling on your original, unrecalibrated data. 
#Then take the SNPs that you have the highest confidence in and use that set as the database of known SNPs
# by feeding it as a VCF file to the base quality score recalibrator. 
#Finally, do a real round of SNP calling with the recalibrated data. 
#These steps could be repeated several times until convergence. 

# Step0: Index reference and dictionary
#cd /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/

#bwa index /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta
#samtools faidx /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta

#java -jar /apps/picard/2.18.17/picard.jar CreateSequenceDictionary \
#reference=/fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#output=/fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.dict

# Step 1: Get a database of known variants

# We chose 3 individuals with highest coverage to call SNPs for our reference database
# sca0700 (Cov=26x) sca0886 (cov =27x) sca1038 (cov =27x)

#cd /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step3/ref_variants/

#samples1=( sca0700 sca0886 sca1038 )

#for line in "${samples1[@]}"
#do
#	gatk --java-options "-Xmx100g -XX:+UseParallelGC -XX:ParallelGCThreads=20" HaplotypeCaller  \
#	-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#	-I /fs/ess/scratch/PAS1533/smathur/wgr/final_bam/old/${line}.Scat.HiCv1.recal.2.bam \
#	-ERC GVCF --native-pair-hmm-threads 20 \
#	-O ${line}.ref.variants.g.vcf 
#done


#gatk --java-options "-Xmx100g -XX:+UseParallelGC -XX:ParallelGCThreads=20" CombineGVCFs \
#-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#--variant sca0700.Scat.HiCv1.raw.variants.g.vcf \
#--variant sca0886.Scat.HiCv1.raw.variants.g.vcf \
#--variant sca1038.Scat.HiCv1.raw.variants.g.vcf \
#-O ref.cohort.g.vcf

#gatk --java-options "-Xmx100g" GenotypeGVCFs \
#-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V ref.cohort.g.vcf \
#-O ref.cohort.vcf 

#gatk SelectVariants \
#-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V ref.cohort.vcf \
#--select-type-to-include SNP \
#-O ref.cohort.snps.vcf
 
#filter variants

#gatk --java-options "-Xmx100g" VariantFiltration \
#-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V ref.cohort.snps.vcf \
#-filter "QD < 2.0" --filter-name "QD2" \
#-filter "ExcessHet > 3.0" --filter-name "ExHet3" \
#-filter "MQ < 40.0" --filter-name "MQ30" \
#-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#-filter "SOR > 3.0" --filter-name "SOR3" \
#-filter "FS > 60.0" --filter-name "FS60" \
#-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#--missing-values-evaluate-as-failing true \
#-select 'vc.isNotFiltered()' \
#-O ref.cohort.snps.filter.vcf

#gatk SelectVariants \
#-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-V ref.cohort.snps.filter.vcf \
#-select 'vc.isNotFiltered()' \
#-O ref.cohort.snps.vcf

# Step 2: Apply BQSR

cd /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step3/

while read -a line
do
	mkdir ${line}
	mkdir ${line}/log/
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/sampleNames.txt


while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}_bqsr
#SBATCH -e ${line}_bqsr
#SBATCH -o ${line}_bqsr

cd $SLURM_SUBMIT_DIR
module load gatk


cd /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step3/${line}/

gatk --java-options \"-Xmx80g\" BaseRecalibrator \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
--known-sites /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step3/ref_variants/ref.cohort.snps.vcf \
-I /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step2/${line}.scat.HiCv1.fixmate.bam \
-O ${line}.recal.1.table &> log/${line}.recal.1.log.txt

gatk --java-options \"-Xmx80g\" ApplyBQSR \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-I /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step2/${line}.scat.HiCv1.fixmate.bam \
--bqsr-recal-file ${line}.recal.1.table \
-O ${line}.Scat.HiCv1.recal.1.bam &> log/${line}.1.ApplyBQSR.log.txt


# ROUND 2 # 

gatk --java-options \"-Xmx80g\" BaseRecalibrator \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
--known-sites /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step3/ref_variants/ref.cohort.snps.vcf \
-I ${line}.Scat.HiCv1.recal.1.bam \
-O ${line}.recal.2.table &> log/${line}.recal.2.log.txt

gatk --java-options \"-Xmx80g\" ApplyBQSR \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-I ${line}.Scat.HiCv1.recal.1.bam \
--bqsr-recal-file ${line}.recal.2.table \
-O ${line}.Scat.HiCv1.recal.2.bam &> log/${line}.2.ApplyBQSR.log.txt" > /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${line}/${line}_baserecal.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/sampleNames.txt













