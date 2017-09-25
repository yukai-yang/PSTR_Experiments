#!/bin/bash
#
#SBATCH --job-name=ECV02_2
#SBATCH --output=output/EXPECV02_2.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1000

Rscript EXPECV02_2.R
