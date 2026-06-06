#!/bin/bash

#SBATCH --job-name=homermmotif
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=homermotif_%j.out
#SBATCH --error=homermotif_%j.err
#SBATCH --partition=compute

# Load system modules
module load anaconda3

# Move to working directory
cd .../Cut_and_Run/

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate homer_env

# Run HOMER motif finding
findMotifsGenome.pl peaks.bed danRer10 homer_motifs/ -size 200 -mask

