#!/bin/sh -e

#SBATCH -p ei-medium #queue
#SBATCH -c 4 #number of cores
#SBATCH --mem 30G #memory pool for all cores
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address

samplename=$1
source trim_galore-0.5.0

trim_galore --paired ${samplename}_R1.fastq.gz ${samplename}_R2.fastq.gz
