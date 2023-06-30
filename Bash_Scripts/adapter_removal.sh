#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=adapter_removal
#SBATCH -e adapter_removal
#SBATCH -o adapter_removal

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 01/26/21                  Last Modified: 06/22/21 ###
###########################################################################
###########################################################################
###                     adapter_removal.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR


#run for each sample

while read -a sample
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=${sample}.adptrem
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
#module load fastqc/0.11.8
module load trimmomatic

# Step1: Adapter removal 

java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
PE /fs/ess/scratch/PAS1533/smathur/wgr/data/raw/${sample}*R1*fastq /fs/ess/scratch/PAS1533/smathur/wgr/data/raw/${sample}*R2*fastq \
/fs/ess/scratch/PAS1533/smathur/wgr/data/trimmed/${sample}.R1.paired /fs/ess/scratch/PAS1533/smathur/wgr/data/trimmed/${sample}.R1.unpaired \
/fs/ess/scratch/PAS1533/smathur/wgr/data/trimmed/${sample}.R2.paired /fs/ess/scratch/PAS1533/smathur/wgr/data/trimmed/${sample}.R2.unpaired \
LEADING:20 TRAILING:20 MINLEN:30 -threads 20 \
ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:40:10


# Step2: FASTQC on raw reads

#cd /fs/scratch/PAS1533/smathur/wgr/read_preprocess/fastqc/trimmed/
#mkdir ${sample}

#fastqc -o ${line}/ \
#/fs/scratch/PAS1533/smathur/wgr/read_preprocess/trimmomatic/trimmed/paired/${sample}*R1*fastq.paired \
#/fs/scratch/PAS1533/smathur/wgr/read_preprocess/trimmomatic/trimmed/paired/${sample}*R2*fastq.paired

#cd /fs/scratch/PAS1533/smathur/wgr/read_preprocess/fastqc/raw/${sample}
#unzip *zip" \
> /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${sample}/${line}.adptrem.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/sampleNames.txt
