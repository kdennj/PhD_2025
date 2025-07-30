#!/bin/bash
#SBATCH -p ei-medium #queue
#SBATCH -c 1 -N 1 #number of cores and nodes
#SBATCH --mem 40G #memory pool for all cores
#SBATCH --array=0-143%15
#SBATCH -J samtools

source samtools-1.7

cd /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/chromosome_coverage_A_M

mapfile -t myID < /KDJ_Scratch_project/KDJ_beans/samples_uniq_all.txt

samtools view -h -b -F 256 ${myID[${SLURM_ARRAY_TASK_ID}]}_pp.bam > primary_aligned/${myID[${SLURM_ARRAY_TASK_ID}]}_primary.bam

