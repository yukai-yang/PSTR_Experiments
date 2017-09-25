#!/bin/bash
#
#SBATCH --job-name=WB106
#SBATCH --output=output/WB106.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWBPOWER106.R

