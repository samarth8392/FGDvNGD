#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH --job-name=aligment
#SBATCH -e aligment_%j
#SBATCH -o aligment_%j

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 01/27/21                  Last Modified: 06/23/21 ###
###########################################################################
###########################################################################
###                     aligment.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR
#module load bwa
#module load samtools/1.8
#module load picard/2.18.17
# step0: Index reference and dictionary

#bwa index /fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta
#samtools faidx /fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta

#java -jar /apps/picard/2.18.17/picard.jar CreateSequenceDictionary \
#reference=/fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#output=/fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.dict


# run for each sample


while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}.align
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load bwa #bwa/0.7.17
module load picard/2.18.17
module load gatk/4.1.2.0
module load samtools/1.8

# STEP_1: Align to Scat genome (HiC_v1.1)

cd /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step1/

bwa mem -t 50 -M -R \"@RG\tID:group1\tSM:${line}\tPL:illumina\tLB:lib1\tPU:unit1\" \
/fs/ess/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
/fs/ess/scratch/PAS1533/smathur/wgr/data/old/${line}.R1.filtered.fastq \
/fs/ess/scratch/PAS1533/smathur/wgr/data/old/${line}.R2.filtered.fastq \
> ${line}.scat.HiCv1.sam

java -jar /apps/picard/2.18.17/picard.jar ValidateSamFile \
I=${line}.scat.HiCv1.sam MODE=SUMMARY O=${line}.scat.HiCv1.sam.txt

java -jar /apps/picard/2.18.17/picard.jar SortSam \
INPUT=${line}.scat.HiCv1.sam OUTPUT=sorted_${line}.scat.HiCv1.bam SORT_ORDER=coordinate

java -jar /apps/picard/2.18.17/picard.jar MarkDuplicates \
INPUT=sorted_${line}.scat.HiCv1.bam OUTPUT=dedup_${line}.scat.HiCv1.bam METRICS_FILE=metrics_${line}.scat.HiCv1.bam.txt

java -jar /apps/picard/2.18.17/picard.jar BuildBamIndex \
INPUT=dedup_${line}.scat.HiCv1.bam

# STEP_2: Fix mate pair info in BAM

cd /fs/ess/scratch/PAS1533/smathur/wgr/alignment/step2/

java -jar /apps/picard/2.18.17/picard.jar FixMateInformation \
INPUT=../step1/dedup_${line}.scat.HiCv1.bam \
OUTPUT=${line}.scat.HiCv1.fixmate.bam \
SO=coordinate \
CREATE_INDEX=true" \
> /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${line}/${line}.align.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/sampleNames.txt

# Get mapping stats

while read -a line
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH --job-name=${line}.align
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge

module load samtools/1.8

# STEP_3: Get mapping stats

cd /fs/ess/scratch/PAS1533/smathur/wgr/alignment/stats/

#samtools depth -a ../step2/round2/${line}.scat.HiCv1.fixmate.bam \
#| awk '{c++;s+=\$3}END{print s/c}' \
#> ${line}.scat.HiCv1.meandepth.txt

#samtools depth -a ../step2/round2/${line}.scat.HiCv1.fixmate.bam \
#| awk '{c++; if(\$3>0) total+=1}END{print (total/c)*100}' \
#> ${line}.scat.HiCv1.1xbreadth.txt

#samtools depth -a ../step2/round2/${line}.scat.HiCv1.fixmate.bam \
#| awk '{c++; if(\$3>=5) total+=1}END{print (total/c)*100}' \
#> ${line}.scat.HiCv1.5xbreadth.txt

samtools depth -a /fs/ess/scratch/PAS1533/smathur/wgr/final_bam/old/${line}.Scat.HiCv1.recal.2.bam \
| awk '{c++; if(\$3>=10) total+=1}END{print (total/c)*100}' \
> ${line}.scat.HiCv1.10xbreadth.txt

#samtools flagstat ../step2/round2/${line}.scat.HiCv1.fixmate.bam \
#> ${line}.scat.HiCv1.mapstats.txt

#samtools stats ../step2/round2/${line}.scat.HiCv1.fixmate.bam \
#> ${line}.scat.HiCv1.stats.txt" \
> /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${line}/${line}.mapstats.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/sampleNames.txt
