#!/bin/bash
#SBATCH --mail-type=END, FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -p ei-long
#SBATCH -c 4
#SBATCH --mem 100G
#SBATCH -J bowtie

cd /KDJ_beans/

sample=$1

source bowtie-1.1.2
source samtools-1.9
source perl-5.22.1

DIR=/KDJ_beans/trim_reads
OUTDIR=/KDJ_beans/bowtie_meso_andean


bowtie -n -q -m 1 -X 450 --fr --al -x ref_genomes/Andean_Meso_ref_genome/Andean_Meso_ref_bowtie \
-1 <(zcat ${DIR}/${sample}_L001_R1_val_1.fq.gz ${DIR}/${sample}_L002_R1_val_1.fq.gz) \
-2 <(zcat ${DIR}/${sample}_L001_R2_val_2.fq.gz ${DIR}/${sample}_L002_R2_val_2.fq.gz) --sam ${sample}_bowtie.sam
#| samtools view -@ 8 -b -S -h - | samtools sort -@ 8 -T ${OUTDIR}/${sample}_tmp -o ${OUTDIR}/${sample}.sorted.bam
