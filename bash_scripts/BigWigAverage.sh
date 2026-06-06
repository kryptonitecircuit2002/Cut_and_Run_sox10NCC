#!/bin/bash

#SBATCH --job-name=bw_avg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=bw_avg_%j.out
#SBATCH --error=bw_avg_%j.err
#SBATCH --partition=compute

#system modules
module load anaconda3

cd .../Cut_and_Run/

#Initialise conda
eval "$(conda shell.bash hook)"

#activate conda environment
conda activate cut_and_run

bigwigAverage -b /.../file1.bw /.../file2.bw /.../file3.bw -o file_norm.bw


