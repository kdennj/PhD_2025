#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -J picard
#SBATCH -c 1 --mem 100G 
#SBATCH --mail-user=denning@nbi.ac.uk
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --array=1-145
#SBATCH -t 1-00:00

source jre-8u144
source picardtools-2.1.1

DIR=/hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bwa_Andean/picardtools

mapfile -t file < samples_uniq_all.txt

java -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar \
MarkDuplicates INPUT=${file[${SLURM_ARRAY_TASK_ID}]}.sorted.bam \
OUTPUT=${DIR}/${file[${SLURM_ARRAY_TASK_ID}]}_A.mdup.bam \
M=${DIR}/${file[${SLURM_ARRAY_TASK_ID}]}_A.mdup.txt
