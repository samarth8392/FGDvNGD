#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=72:00:00
#SBATCH --account=PAS1533
#SBATCH --mail-user=mason.501@osu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=Aug_Train2_Scate_HiC



cd $SLURM_SUBMIT_DIR

source activate maker

gff3_merge -n -s -d /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run2/Scate_HiC_rnd2.maker.output/Scate_HiC_rnd2_master_datastore_index.log > /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run2/Scate_HiC_rnd2.maker.output/Scate_HiC_rnd2.all.noseq.gff



conda deactivate
source activate chimerakiller_env

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' Scate_HiC_rnd2.all.noseq.gff | awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | bedtools getfasta -fi /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Scatenatus_HiC_v1.1.fasta -bed - -fo /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run2/Scate_HiC_rnd2.all.maker.transcripts1000.fasta

echo "finished"
