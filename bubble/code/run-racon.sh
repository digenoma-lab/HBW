#!/bin/bash
#SBATCH --job-name="Racon"
#SBATCH --partition=general
#SBATCH -c 44
#SBATCH --mem-per-cpu=4200


make -f polishing.mk all CPU=44 

