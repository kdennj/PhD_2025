#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 4
#SBATCH --mem 50G
#SBATCH --array=0-145
#SBATCH --mail-user=denning@nbi.ac.uk
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH -o /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bwa_Meso/slurm.%N.%j.out # STDOUT
#SBATCH -e /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bwa_Meso/slurm.%N.%j.err # STDERR
#SBATCH -J bwamem

source samtools-1.7
source bwa-0.7.13

REFDIR=/ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/ref_genomes/Andean_assembly
DIR=/ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/trim_reads
OUTDIR=/ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/bwa_Andean

mapfile -t file < ../samples_uniq_all.txt

bwa mem -M -R "@RG\tID:${file[${SLURM_ARRAY_TASK_ID}]}\tSM:${file[${SLURM_ARRAY_TASK_ID}]}\tLB:${file[${SLURM_ARRAY_TASK_ID}]}\tPL:ILLUMINA\tPU:EI" \
${REFDIR}/Pvulgaris_442_v2.0.fa \
<(zcat ${DIR}/${file[${SLURM_ARRAY_TASK_ID}]}_L001_R1_val_1.fq.gz ${DIR}/${file[${SLURM_ARRAY_TASK_ID}]}_L002_R1_val_1.fq.gz) \
<(zcat ${DIR}/${file[${SLURM_ARRAY_TASK_ID}]}_L001_R2_val_2.fq.gz ${DIR}/${file[${SLURM_ARRAY_TASK_ID}]}_L002_R2_val_2.fq.gz) \
| samtools view -@ 8 -b -S -h - | samtools sort -@ 8 -T ${OUTDIR}/${file[${SLURM_ARRAY_TASK_ID}]}_tmp \
-o ${OUTDIR}/${file[${SLURM_ARRAY_TASK_ID}]}.sorted.bam
