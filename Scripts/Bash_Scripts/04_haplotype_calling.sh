#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=variant_calling
#SBATCH -e variant_calling
#SBATCH -o variant_calling

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 02/07/21                  Last Modified: 07/29/21 ###
###########################################################################
###########################################################################
###                     variant_calling.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}_hapCall
#SBATCH -e ${line}_hapCall
#SBATCH -o ${line}_hapCall

cd $SLURM_SUBMIT_DIR
module load gatk

cd /fs/ess/scratch/PAS1533/smathur/wgr/gatk/gvcf/

gatk --java-options \"-Xmx100g -XX:+UseParallelGC -XX:ParallelGCThreads=20\" HaplotypeCaller  \
-R /fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
-I /fs/ess/scratch/PAS1533/smathur/wgr/bamfiles/final/${line}.Scat.HiCv1.recal.2.bam \
-ERC GVCF --native-pair-hmm-threads 20 \
-O ${line}.Scat.HiCv1.raw.variants.g.vcf"  > /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${line}/${line}_HapCall.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/sampleNames.txt
