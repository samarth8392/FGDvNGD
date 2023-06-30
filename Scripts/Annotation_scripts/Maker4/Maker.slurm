#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --account=PAS1533
#SBATCH --mail-user=mason.501@osu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=Scate_HiC_Maker_rnd4


export  REPEATMASKER_DIR=/users/PAS1533/ajmason/miniconda3/envs/maker/share/RepeatMasker
export  REPEATMASKER_LIB_DIR=/users/PAS1533/ajmason/miniconda3/envs/maker/share/RepeatMasker/Libraries
export  REPEATMASKER_MATRICES_DIR=/users/PAS1533/ajmason/miniconda3/envs/maker/share/RepeatMasker/Matrices


cd $SLURM_SUBMIT_DIR


source activate maker

mpiexec -n 6 maker -fix_nucleotides -base Scate_HiC_rnd4 maker_bopts.ctl maker_exe.ctl maker_opts.ctl

echo "finished"
