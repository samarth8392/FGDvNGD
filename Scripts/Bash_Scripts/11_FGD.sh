#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 150:00:00
#SBATCH --job-name=FGD
#SBATCH -e %x
#SBATCH -o %x

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 09/13/21                  Last Modified: 05/12/23 ###
###########################################################################
###########################################################################
###                     fgd.sh                        					###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load vcftools
module load samtools
module load gcc-compatibility
module load htslib

cd /fs/scratch/PAS1533/smathur/wgr/fgd/


################# 
# MUTATION LOAD #
################# 

############################################################ 

# STEP 1a: Get population frequency and individual genotypes
# (Using all samples)

while read -a pop
do
	for sites in lof.lower10 upper10
	do
		vcftools --vcf /fs/scratch/PAS1533/smathur/wgr/vcf/bypop/${pop}/onlyScat.only_${pop}.wholegenome.SNPs.recode.vcf \
		--positions /fs/scratch/PAS1533/smathur/wgr/sites/${sites}.sites.txt \
		--freq \
		--out freq/${pop}.${sites}

		sed -i '1d' freq/${pop}.${sites}.frq
		less freq/${pop}.${sites}.frq | cut -f6 | cut -f2 -d ":" > freq/${pop}.${sites}.DAF
		paste freq/${pop}.${sites}.frq freq/${pop}.${sites}.DAF > freq/${pop}.${sites}.freq

		rm freq/*.frq
		rm freq/*.DAF

		vcftools --vcf /fs/scratch/PAS1533/smathur/wgr/vcf/bypop/${pop}/onlyScat.only_${pop}.wholegenome.SNPs.recode.vcf \
		--positions /fs/scratch/PAS1533/smathur/wgr/sites/${sites}.sites.txt \
		--extract-FORMAT-info GT \
		--out geno/${pop}.${sites}
	done
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames_noGRLL.txt


# STEP 1b: Bootstrapping 100 times using only 5 samples per population

cd /fs/scratch/PAS1533/smathur/wgr/fgd/boot/

while read -a pop
do
	cd /fs/scratch/PAS1533/smathur/wgr/fgd/boot/
	mkdir ${pop}
	for r in {1..100}
	do
		shuf -n 5 /fs/scratch/PAS1533/smathur/wgr/lists/${pop}.sampleList > ${pop}/${pop}.${r}.samples
		for sites in lof.lower10 upper10
		do
			vcftools --vcf /fs/scratch/PAS1533/smathur/wgr/vcf/bypop/${pop}/onlyScat.only_${pop}.wholegenome.SNPs.recode.vcf \
			--positions /fs/scratch/PAS1533/smathur/wgr/sites/${sites}.sites.txt \
			--freq \
			--keep ${pop}/${pop}.${r}.samples \
			--out ${pop}/${pop}.run${r}.${sites}

			sed -i '1d' ${pop}/${pop}.run${r}.${sites}.frq
			less ${pop}/${pop}.run${r}.${sites}.frq | cut -f6 | cut -f2 -d ":" > ${pop}/${pop}.run${r}.${sites}.DAF
			paste ${pop}/${pop}.run${r}.${sites}.frq ${pop}/${pop}.run${r}.${sites}.DAF > ${pop}/${pop}.run${r}.${sites}.freq

			rm ${pop}/*.frq
			rm ${pop}/*.DAF
		done
	done
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames_noGRLL.txt


# STEP 2: Estimate mutation load using Rscript geneticLoad.R

############################################################ 


###################### 
# ADAPTIVE DIVERSITY #
###################### 

# Upper10% Pi Using ANGSD

# index sites

cd /fs/ess/scratch/PAS1533/smathur/wgr/sites/
~/softwares/angsd/angsd sites index upper10.sites.txt


# index bam files
cd /fs/ess/scratch/PAS1533/smathur/wgr/geneDes/bams/inGene/
for sample in $(ls *bam| cut -f1 -d ".")
do
	samtools index -b ${sample}.inGene.bam
done


while read -a pop
do
	echo "#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 150:00:00
#SBATCH --job-name=${pop}.up10
#SBATCH -e %x
#SBATCH -o %x

cd $SLURM_SUBMIT_DIR
module purge
module load samtools
module load gcc-compatibility
module load htslib

cd /fs/ess/scratch/PAS1533/smathur/wgr/pi/upper10/

# Get 1D SFS

#mkdir ${pop}

#~/softwares/angsd/angsd -P 128 \
#-dosaf 1 -GL 1 \
#-anc /fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-bam /fs/scratch/PAS1533/smathur/wgr/lists/bamlist/${pop}.inGene.bamlist \
#-rf /fs/ess/scratch/PAS1533/smathur/wgr/sites/upper10_justCDS.regions.txt \
#-out ${pop}/${pop}.upper10_justCDS 

#~/softwares/angsd/misc/realSFS -P 128 -fold 1 \
#${pop}/${pop}.upper10_justCDS.saf.idx -maxIter 100 \
#> ${pop}/${pop}.upper10_justCDS.sfs

# Get thetas

#~/softwares/angsd/misc/realSFS saf2theta \
#${pop}/${pop}.upper10_justCDS.saf.idx \
#-sfs ${pop}/${pop}.upper10_justCDS.sfs -fold 1 \
#-outname ${pop}/${pop}.upper10_justCDS

#~/softwares/angsd/misc/thetaStat do_stat \
#${pop}/${pop}.upper10_justCDS.thetas.idx \
#-outnames ${pop}/${pop}.upper10_justCDS

#~/softwares/angsd/misc/thetaStat print \
#${pop}/${pop}.upper10_justCDS.thetas.idx \
#> ${pop}/${pop}.upper10_justCDS_persite.txt

# Dump Counts
~/softwares/angsd/angsd -P 128 \
-doCounts 1 -dumpCounts 2 \
-rf /fs/ess/scratch/PAS1533/smathur/wgr/sites/upper10_justCDS.regions.txt \
-bam /fs/scratch/PAS1533/smathur/wgr/lists/bamlist/${pop}.inGene.bamlist \
-out ${pop}/${pop}.upper10_justCDS 


# With sites

#~/softwares/angsd/angsd -P 128 \
#-dosaf 1 -GL 1 \
#-anc /fs/scratch/PAS1533/smathur/wgr/ref/scat/Scatenatus_HiC_v1.1.fasta \
#-bam /fs/scratch/PAS1533/smathur/wgr/lists/bamlist/${pop}.inGene.bamlist \
#-sites /fs/ess/scratch/PAS1533/smathur/wgr/sites/upper10.sites.txt \
#-out ${pop}/${pop}.upper10_sites 

#~/softwares/angsd/misc/realSFS -P 128 -fold 1 \
#${pop}/${pop}.upper10_sites.saf.idx -maxIter 100 \
#> ${pop}/${pop}.upper10_sites.sfs

# Get thetas

#~/softwares/angsd/misc/realSFS saf2theta \
#${pop}/${pop}.upper10_sites.saf.idx \
#-sfs ${pop}/${pop}.upper10_sites.sfs -fold 1 \
#-outname ${pop}/${pop}.upper10_sites

#~/softwares/angsd/misc/thetaStat do_stat \
#${pop}/${pop}.upper10_sites.thetas.idx \
#-outnames ${pop}/${pop}.upper10_sites

#~/softwares/angsd/misc/thetaStat print \
#${pop}/${pop}.upper10_sites.thetas.idx \
#> ${pop}/${pop}.upper10_persite.txt

# Dump Counts
~/softwares/angsd/angsd -P 128 \
-doCounts 1 -dumpCounts 2 \
-sites /fs/ess/scratch/PAS1533/smathur/wgr/sites/upper10.sites.txt \
-bam /fs/scratch/PAS1533/smathur/wgr/lists/bamlist/${pop}.inGene.bamlist \
-out ${pop}/${pop}.upper10_sites" \
> /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/perPop/${pop}.upper10.theta.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames_noGRLL.txt

# Submit jobs

while read -a pop
do
	cd /fs/ess/scratch/PAS1533/smathur/wgr/errors/perPop/
	sbatch /fs/ess/scratch/PAS1533/smathur/wgr/jobcodes/perPop/${pop}.upper10.theta.sh
done < /fs/ess/scratch/PAS1533/smathur/wgr/lists/popnames_noGRLL.txt