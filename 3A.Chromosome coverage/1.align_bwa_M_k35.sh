#!/bin/bash -e
#SBATCH -p ei-long
#SBATCH -N 1 # number of nodes
#SBATCH -J bwa-mem
#SBATCH -c 1 # number of tasks
#SBATCH --mem 30G # memory pool for all cores

#- @ 8 would make run faster in the future
mysamplename=$1

DIR=/KDJ_beans_project/trim_reads
OUTDIR=/KDJ_beans_project/bwa_meso_andean

source bwa-0.7.17
source samtools-1.7

bwa mem -M -k 35 -t 1 -R "@RG\tID:${mysamplename}\tSM:${mysamplename}\tPL:ILLUMINA\tPU:EI" \
/ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans_project/ref_genomes/Andean_Meso_ref.fa \
<(zcat ${DIR}/${mysamplename}_L001_R1_val_1.fq.gz ${DIR}/${mysamplename}_L002_R1_val_1.fq.gz) \
<(zcat ${DIR}/${mysamplename}_L001_R2_val_2.fq.gz ${DIR}/${mysamplename}_L002_R2_val_2.fq.gz) \
| samtools view -@ 8 -b -S -h - | samtools sort -@ 8 -T ${OUTDIR}/${mysamplename}_tmp -o ${OUTDIR}/${mysamplename}.sorted.bam

