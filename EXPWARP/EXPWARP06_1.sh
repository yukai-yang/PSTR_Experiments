#!/bin/bash
#
#SBATCH --job-name=WARP06_1
#SBATCH --output=output/EXPWARP06_1.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP06_1.R
