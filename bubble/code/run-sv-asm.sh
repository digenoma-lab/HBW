#!/bin/bash
#SBATCH --job-name="SV-ASM"
#SBATCH --partition=slims
#SBATCH --share
#SBATCH -c 20
#SBATCH --mem-per-cpu=2000

ml Anaconda3/5.3.0 squashfs-tools/4.4
CONTAINER=/home/adigenova/HSCAFF/CANCER/BENCHMARK/HG0002/ont_callers_v1.0.sif
REF=/home/adigenova/HSCAFF/CANCER/BENCHMARK/HG0002/ref/human_hs37d5.fasta
#singularity exec ${CONTAINER} make -f sv-asm.mk REF=${REF} CPU=20 PREF=WGD_HG002 HAP1=wengan.hg002.hap1.fa HAP2=wengan.hg002.hap2.fa  all 
#singularity exec ${CONTAINER} make -f sv-asm.mk REF=${REF} CPU=20 PREF=HF_HG002 HAP1=NA24385.HiFi.hifiasm-0.12.hap1.fa.gz HAP2=NA24385.HiFi.hifiasm-0.12.hap2.fa.gz  all 
singularity exec ${CONTAINER} make -f sv-asm.mk REF=${REF} CPU=20 PREF=WGD_HG002-R2 HAP1=asm1.r3.fa HAP2=asm1.r3.hap2.fa all
singularity exec ${CONTAINER} make -f sv-asm-def.mk REF=${REF} CPU=20 PREF=WGD_HG002-R2-DEF HAP1=asm1.r3.fa HAP2=asm1.r3.hap2.fa all
#make -f  asm1.mk all 
