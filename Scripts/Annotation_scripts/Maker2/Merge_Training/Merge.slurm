#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=168:00:00
#SBATCH --account=PAS1533
#SBATCH --mail-user=mason.501@osu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=Merge_Maker2_Scate_HiC


cd $SLURM_SUBMIT_DIR


export  REPEATMASKER_DIR=/users/PAS1533/ajmason/miniconda3/envs/maker/share/RepeatMasker
export  REPEATMASKER_LIB_DIR=/users/PAS1533/ajmason/miniconda3/envs/maker/share/RepeatMasker/Libraries
export  REPEATMASKER_MATRICES_DIR=/users/PAS1533/ajmason/miniconda3/envs/maker/share/RepeatMasker/Matrices


source activate maker

fasta_merge -d /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run2/Scate_HiC_rnd2.maker.output/Scate_HiC_rnd2_master_datastore_index.log

gff3_merge -d /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run2/Scate_HiC_rnd2.maker.output/Scate_HiC_rnd2_master_datastore_index.log

echo "finished"
