#!/bin/bash
#
#SBATCH --job-name=WB201
#SBATCH --output=output/WB201.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=70:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPWBPOWER201.R

