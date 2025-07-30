#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH -c 1 --mem 100G 
#SBATCH -J mkdup # number of tasks
#SBATCH --mem 60G # memory pool for all cor 

mySAMN=$1
source jre-8u144
source picardtools-2.1.1

cd /scratch/KDJ_beans/bowtie2_meso_andean


java -jar /tgac/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar \
MarkDuplicates INPUT=sorted_bams/${mySAMN}.sorted.bam \
OUTPUT= pp/${mySAMN}.mdup.bam \
M= pp/${mySAMN}.mdup.txt REMOVE_DUPLICATES=true
