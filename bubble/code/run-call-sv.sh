#!/bin/bash
#SBATCH --job-name="long"
#SBATCH --partition=largemem
#SBATCH -c 44
#SBATCH --mem-per-cpu=4000


#make all CPU=44 
make -f call_sv_asm.mk CPU=44 PREF=wengan.r3 REF=asm1.r3.fa all

