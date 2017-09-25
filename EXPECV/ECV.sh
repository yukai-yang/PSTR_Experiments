#!/bin/bash
#
#SBATCH --job-name=ECV
#SBATCH --output=output/ECV.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=1000

Rscript ECV.R 
