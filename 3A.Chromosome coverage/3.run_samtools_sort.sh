#!/bin/bash
#SBATCH -p ei-medium #queue
#SBATCH -c 4 #number of cores
#SBATCH -N 1 # number of nodes
#SBATCH --mem 30G #memory pool for all cores
#SBATCH -J samtools
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -o /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bowtie2_meso_andean/sorted_bams/slurm.%N.%j.out # STDOUT

sample=$1

source samtools-1.7


###Script for Andean_Meso bowtie2
# cd /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bowtie2_meso_andean
# samtools view -@ 8 -b -S -h ${sample}_bowtie2.bam | samtools sort -@ 8 -T sorted_bams/${sample}_tmp -o sorted_bams/${sample}.sorted.bam


#Script for Andean bwa
#samtools sort -m 5G -O bam -T ${sample}_tmp -o ${sample}.sorted.bam /ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans_project/bwa_Andean_M/${sample}.bam
#samtools view -b -S -h ./${sample} | samtools sort -T ${sample}_temporal > ${sample}.sorted.bam
#| samtools view -@ 8 -b -S -h - | samtools sort -@ 8 -T ${OUTDIR}/${sample}_tmp -o ${OUTDIR}/${sample}.sorted.bam
