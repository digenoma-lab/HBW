#!/bin/bash
#SBATCH --job-name="long"
#SBATCH --partition=general
#SBATCH -c 44
#SBATCH --mem-per-cpu=4200


make all CPU=44 

