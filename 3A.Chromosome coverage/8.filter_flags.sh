#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -c 1 # number of tasks
#SBATCH --mem 150G # memory pool for all cores

myID=$1

source samtools-1.7

cd /ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/bowtie2_meso_andean/pp

samtools view -h -f 0x2 ${myID}.mdup.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ../unique/${myID}_unique.bam