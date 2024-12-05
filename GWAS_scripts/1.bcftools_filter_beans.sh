#!/bin/bash -e
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -J bcftools
#SBATCH -p ei-medium
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH --mem 20G # memory pool for all cores


#Filter for depth
#Filter for bi-allelic SNPs / maximum alleles 2

source bcftools-1.12

myID=$1

bcftools filter -e "FMT/DP<5" -S "." ${myID}.vcf | \
bcftools filter -i 'QUAL>30' | \
bcftools view -i 'F_MISSING<0.1' | \ #Max Missing
bcftools view -m2 -M2 -O v > /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/bwa_Andean_M/gatk_filter_GWAS/${myID}_filtered.vcf.biallelic_MAF01_thin5bp.vcf.gz
