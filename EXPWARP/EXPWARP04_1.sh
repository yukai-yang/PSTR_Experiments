#!/bin/bash
#
#SBATCH --job-name=WARP04_1
#SBATCH --output=output/EXPWARP04_1.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP04_1.R
