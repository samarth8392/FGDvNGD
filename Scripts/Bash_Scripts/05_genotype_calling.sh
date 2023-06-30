#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 150:00:00
#SBATCH --job-name=genotypeGVCF_all
#SBATCH --mem=2800G
#SBATCH -e genotypeGVCF_all
#SBATCH -o genotypeGVCF_all

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 02/07/21                  Last Modified: 08/04/21 ###
###########################################################################
###########################################################################
###                     genotypeGVCF.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load gatk #gatk/4.1.2.0
module load vcftools

## Step 1: CombineGVCFs 

cd /fs/ess/scratch/PAS1533/smathur/wgr/gatk/gvcf/

gatk --java-options "-Xmx2800g -XX:+UseParallelGC -XX:ParallelGCThreads=48" CombineGVCFs  \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
--variant sca0142.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0143.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0186.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0226.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0236.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0254.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0260.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0673.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0678.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0689.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0693.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0694.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0700.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0707.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0708.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0709.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0713.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0724.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0729.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0732.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0734.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0742.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0744.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0750.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0801.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0802.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0807.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0808.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0809.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0810.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0812.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0863.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0884.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0885.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0886.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0887.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0888.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0893.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0924.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0925.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0926.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0931.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0932.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0933.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0939.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0940.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0965.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0966.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0967.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0968.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0978.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0979.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0980.Scat.HiCv1.raw.variants.g.vcf \
--variant sca0982.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1018.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1019.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1020.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1023.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1026.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1027.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1028.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1029.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1030.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1035.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1036.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1037.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1038.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1050.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1055.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1057.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1058.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1078.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1223.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1224.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1227.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1228.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1229.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1230.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1231.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1259.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1266.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1272.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1372.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1379.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1381.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1395.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1453.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca1515.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1529.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1530.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1531.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1532.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1536.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1541.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1542.Scat.HiCv1.raw.variants.g.vcf \
--variant sca1543.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca461.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca465.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca472.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca529.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca531.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca549.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca600.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca604.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca636.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca651.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca901.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca903.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca905.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca917.Scat.HiCv1.raw.variants.g.vcf \
--variant Sca919.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0005.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0096.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0097.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0098.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0105.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0111.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0114.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0119.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0120.Scat.HiCv1.raw.variants.g.vcf \
--variant ste0121.Scat.HiCv1.raw.variants.g.vcf \
-O all121.final.10x.g.vcf.gz



# Step2: GenotypeGVCF

cd /fs/ess/scratch/PAS1533/smathur/wgr/gatk/variants/

gatk --java-options "-Xmx2800g -XX:+UseParallelGC -XX:ParallelGCThreads=48" GenotypeGVCFs \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-V /fs/ess/scratch/PAS1533/smathur/wgr/gatk/gvcf/all121.final.10x.g.vcf.gz \
-O all121.final.10x.raw.variants.vcf.gz


# Step3: Select SNPs

gatk --java-options "-Xmx2800g -XX:+UseParallelGC -XX:ParallelGCThreads=48" SelectVariants \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-V all121.final.10x.raw.variants.vcf.gz \
--select-type-to-include SNP \
-O all121.final.10x.raw.SNPs.vcf.gz 

#Filter SNPs
# From grossen et al 2020 QD < 2.0, FS > 40.0, SOR > 5.0, MQ< 20.0, −3.0 > MQRandkSum >3.0, 
# −3.0 > ReadPosRankSum >3.0 and AN < 62 (80% of all Alpine ibex individuals)

# using AN < 193 (80% of all samples (2X121 alleles))


gatk --java-options "-Xmx2800g" VariantFiltration \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-V all121.final.10x.raw.SNPs.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "MQ < 20.0" --filter-name "MQ20" \
-filter "MQRankSum < -3.0" --filter-name "MQRankSum-3.0" \
-filter "MQRankSum > 3.0" --filter-name "MQRankSum3.0" \
-filter "ReadPosRankSum < -3.0" --filter-name "ReadPosRankSum-3.0" \
-filter "ReadPosRankSum > 3.0" --filter-name "ReadPosRankSum3.0" \
-filter "SOR > 5.0" --filter-name "SOR5" \
-filter "FS > 40.0" --filter-name "FS40" \
-filter "AN < 193.0" --filter-name "AN193" \
-filter "AF < 0.05" --filter-name "MAF0.05" \
-O all121.final.10x.filterflag.SNPs.vcf.gz

gatk SelectVariants \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-V all121.final.10x.filterflag.SNPs.vcf.gz \
-select 'vc.isNotFiltered()' \
-O all121.final.10x.filtered.SNPs.vcf.gz

# keep only biallelic sites

vcftools --gzvcf all121.final.10x.filtered.SNPs.vcf.gz \
--recode --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 \
--out all121.final.10x.final.SNPs.vcf.gz


