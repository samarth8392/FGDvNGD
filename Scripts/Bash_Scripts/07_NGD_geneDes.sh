#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=geneDes_pi
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/13/21                  Last Modified: 04/24/23 ###
###########################################################################
###########################################################################
###                     geneDes_pi.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load samtools

# Step0 : Reindex fasta

cd /fs/scratch/PAS1533/smathur/wgr/ref/scat/
samtools faidx Scatenatus_HiC_v1.1.fasta

# Step1: Remove all genic regions from BAM file 
# geneRegions.bed and 500kb.updown.genes.bed are provided in data folder

cd /fs/ess/scratch/PAS1533/smathur/wgr/bamfiles/

for sample in $(ls *bam| cut -f1 -d ".")
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 150:00:00
#SBATCH --job-name=${sample}.bamFilt
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load samtools
module load htslib

cd /fs/ess/scratch/PAS1533/smathur/wgr/geneDes/bams/

samtools view /fs/ess/scratch/PAS1533/smathur/wgr/bamfiles/${sample}.Scat.HiCv1.recal.2.bam \
-b -h -o ${sample}.inGene.bam -U ${sample}.outGene.bam \
-L /fs/ess/scratch/PAS1533/smathur/wgr/geneDes/geneRegions.bed 

samtools view ${sample}.outGene.bam \
-b -h -o ${sample}.inGene500.bam -U ${sample}.outGene500.bam \
-L /fs/ess/scratch/PAS1533/smathur/wgr/geneDes/500kb.updown.genes.bed" \
> /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/perSample/${sample}.bamFilt.sh
done

# Submit jobs
cd /fs/ess/scratch/PAS1533/smathur/wgr/bamfiles/

for sample in $(ls *bam| cut -f1 -d ".")
do
	cd /fs/ess/scratch/PAS1533/smathur/wgr/errors/perSample/
	sbatch  /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/perSample/${sample}.bamFilt.sh
done

# Step2: Run ANGSD for each population separately

while read -a pop
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=${pop}.sfs
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load samtools
module load gcc-compatibility
module load htslib

cd /fs/scratch/PAS1533/smathur/wgr/pi/geneDes/

# Get 1D SFS

#mkdir ${pop}

#~/softwares/angsd/angsd -P 128 \
#-dosaf 1 -GL 1 \
#-anc /fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-bam /fs/scratch/PAS1533/smathur/wgr/lists/bamlist/${pop}.outGene500.bamlist \
#-out ${pop}/${pop}.outGene500 

~/softwares/angsd/misc/realSFS -P 128 -fold 1 \
${pop}/${pop}.outGene500.saf.idx -maxIter 100 \
> ${pop}/${pop}.outGene500.sfs

# Get thetas

~/softwares/angsd/misc/realSFS saf2theta \
${pop}/${pop}.outGene500.saf.idx \
-sfs ${pop}/${pop}.outGene500.sfs -fold 1 \
-outname ${pop}/${pop}.outGene500

~/softwares/angsd/misc/thetaStat do_stat \
${pop}/${pop}.outGene500.thetas.idx \
-win 1000 -step 500 \
-outnames ${pop}/${pop}.outGene500_1_0.5

~/softwares/angsd/misc/thetaStat do_stat \
${pop}/${pop}.outGene500.thetas.idx \
-win 10000 -step 5000 \
-outnames ${pop}/${pop}.outGene500_10_5

~/softwares/angsd/misc/thetaStat do_stat \
${pop}/${pop}.outGene500.thetas.idx \
-outnames ${pop}/${pop}.outGene500

~/softwares/angsd/misc/thetaStat print \
${pop}/${pop}.outGene500.thetas.idx \
> ${pop}/${pop}.outGene500_persite.txt" \
> /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/perPop/${pop}.theta.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames_noGRLL.txt


# Submit jobs
cd /fs/ess/scratch/PAS1533/smathur/wgr/bamfiles/

while read -a pop
do
	cd /fs/ess/scratch/PAS1533/smathur/wgr/errors/perPop/
	sbatch /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/perPop/${pop}.theta.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames_noGRLL.txt

