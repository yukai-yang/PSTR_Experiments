#!/bin/bash
#
#SBATCH --job-name=WARP05_1
#SBATCH --output=output/EXPWARP05_1.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP05_1.R
