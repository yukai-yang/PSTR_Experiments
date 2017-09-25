#!/bin/bash
#
#SBATCH --job-name=R1WB03
#SBATCH --output=output/R1WB03.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWBSIZE03.R
