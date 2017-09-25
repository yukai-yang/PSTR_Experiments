#!/bin/bash
#
#SBATCH --job-name=WARP01
#SBATCH --output=output/EXPWARP01.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP01.R
