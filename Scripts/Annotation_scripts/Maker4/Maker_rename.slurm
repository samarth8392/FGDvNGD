#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --account=PAS1533
#SBATCH --mail-user=mason.501@osu.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=Scate_HiC_rename


cd $SLURM_SUBMIT_DIR

source activate maker
#Create noseq file
gff3_merge -n -s -d /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.maker.output/Scate_HiC_rnd4_master_datastore_index.log > Scate_HiC_rnd4.all.noseq.gff
#Create naming table
maker_map_ids --prefix Scate --justify 5 /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.maker.output/Scate_HiC_rnd4.all.gff > /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.all.maker.name.map
#Replace names in GFF files
map_gff_ids /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.all.maker.name.map /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.maker.output/Scate_HiC_rnd4.all.gff 

map_gff_ids /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.all.maker.name.map /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.all.noseq.gff

#Replace names in FASTA headers
map_fasta_ids /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.all.maker.name.map /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.maker.output/Scate_HiC_rnd4.all.maker.transcripts.fasta

map_fasta_ids /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.all.maker.name.map /fs/scratch/PAS1533/andrewmason/Annotation/Scate_HiC/Maker/Run4/Scate_HiC_rnd4.maker.output/Scate_HiC_rnd4.all.maker.proteins.fasta

echo "finished"

