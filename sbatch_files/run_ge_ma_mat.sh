#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60GB
#SBATCH --job-name="GPU Coloc GE MA mat creation"
#SBATCH --partition=amd

module load python/3.11.1
module load py-pandas/2.0.1
module load py-tqdm/4.64.0

python3 ge_ma_mat.py
