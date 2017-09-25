#!/bin/bash
#
#SBATCH --job-name=WARP07_2
#SBATCH --output=output/EXPWARP07_2.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP07_2.R
