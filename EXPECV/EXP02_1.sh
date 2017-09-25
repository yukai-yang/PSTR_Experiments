#!/bin/bash
#
#SBATCH --job-name=ECV02_1
#SBATCH --output=output/EXPECV02_1.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPECV02_1.R 
