#!/bin/bash -e
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -J bcftools
#SBATCH -p ei-medium
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH --mem 20G # memory pool for all cores

#Filter for depth
#Max Missing
#Keep only SNPs
#Filter for bi-allelic SNPs / maximum alleles 2
#MAF for ref or alternative


source bcftools-1.12

myID=$1

bcftools filter -e "FMT/DP<5" -S "." ${myID}.vcf | \
bcftools filter -i 'QUAL>30' | \
bcftools view -i 'F_MISSING<0.95' |\
bcftools view --types snps | \
bcftools view -m2 -M2 | \
bcftools view -q 0.02:minor -O v > /hpc-home/denning/KDJ_Scratch_project/KDJ_beans_project/bwa_Andean/gatk_filter_QUAL/${myID}_out_filtered_QUAL.vcf