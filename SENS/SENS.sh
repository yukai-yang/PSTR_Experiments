#!/bin/bash
#
#SBATCH --job-name=SENS
#SBATCH --output=output/SENS.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1000

Rscript SENS.R 
