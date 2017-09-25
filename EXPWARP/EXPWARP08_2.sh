#!/bin/bash
#
#SBATCH --job-name=WARP08_2
#SBATCH --output=output/EXPWARP08_2.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP08_2.R 
