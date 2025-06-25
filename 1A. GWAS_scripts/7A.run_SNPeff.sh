#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -J SNPeff
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH --mem 40G # memory pool for all cores

source snpeff-4.3
source jre-8u144

###################
#Step 1 - create database from gff3 and fasta of reference genome 
#Andean genome first
java -Xmx40g -jar /ei/software/testing/snpeff/4.3q/x86_64/snpEff/snpEff.jar build -gff3 -v Andean_ref_genome