#!/bin/bash
#
#SBATCH --job-name=WARP02_2
#SBATCH --output=output/EXPWARP02_2.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWARP02_2.R
