#!/bin/bash
#SBATCH --job-name="Racon"
#SBATCH --partition=largemem
#SBATCH -c 44
#SBATCH --mem-per-cpu=17000


make -f polishing.mk all CPU=44 

