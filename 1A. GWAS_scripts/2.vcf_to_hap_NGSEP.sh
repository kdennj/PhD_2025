#!/bin/bash
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -J NGSEP
#SBATCH -p ei-medium
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH --mem 50G # memory pool for all cores

source jdk-13.0.2
source NGSEP-4.2.1

#Hapmaps to vcfs
java -Xmx40g -jar /ei/software/testing/NGSEP/4.2.1/x86_64/bin/NGSEPcore_4.2.1.jar \
VCFConverter -hapmap -i ${myID}_filtered.vcf.biallelic_MAF01_thin5bp.vcf.gz -o ${myID}_filtered.vcf.biallelic_MAF01_thin5bp.vcf.gz.hapmap.txt_hmp.txt