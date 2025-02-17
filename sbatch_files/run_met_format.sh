#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60GB
#SBATCH --job-name="GPU Coloc format metabolites"
#SBATCH --partition=amd
#SBATCH --mail-user=mihkelje@ut.ee
#SBATCH --mail-type=BEGIN,END,FAIL

module load python/3.11.1
source venv_met/bin/activate
python3 met_format.py
