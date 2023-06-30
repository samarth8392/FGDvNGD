#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --account=PAS1533
#SBATCH --mail-user=mason.501@osu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=BUSCO_Scate_HiC_Augustus_Training


cd $SLURM_SUBMIT_DIR


source activate busco

export BUSCO_CONFIG_FILE="/fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Augustus_Training_2/config.ini"

python /users/PAS1533/ajmason/miniconda3/envs/busco/bin/run_BUSCO.py -r -i /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run2/Scate_HiC_rnd2.all.maker.transcripts1000.fasta  -o Scate_HiC_rnd2_Aug_train -l /users/PAS1533/ajmason/tetrapoda_odb9/ -m genome -c 12 --long -r -sp Sistrurus_catenatus_HiC -z --augustus_parameters='--progress=true'

source deactivate

echo "finished"
