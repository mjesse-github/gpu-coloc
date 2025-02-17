#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="GPU Coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-user=mihkelje@ut.ee
#SBATCH --mail-type=BEGIN,END,FAIL

module load python/3.11.1
source torch_env/bin/activate
python3 coloc.py
