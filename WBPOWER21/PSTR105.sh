#!/bin/bash
#
#SBATCH --job-name=WB105
#SBATCH --output=output/WB105.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWBPOWER105.R

