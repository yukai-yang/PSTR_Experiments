#!/bin/bash
#
#SBATCH --job-name=WB103
#SBATCH --output=output/WB103.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWBPOWER103.R 

