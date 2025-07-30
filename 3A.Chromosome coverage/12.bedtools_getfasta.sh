#!/bin/bash -e
#SBATCH -p ei-short
#SBATCH -J bedtools
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH --mem 20G # memory pool for all cores

cd /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/ref_genomes/

source bedtools-2.26.0
source minimap2-2.22_CBG

bedtools getfasta -fi Andean_assembly/Pvulgaris_Andean_442_v2.0.fa -bed Andean_Meso_ref_genome/AM_window100kb_rep.bed -fo Pvulgaris_Andean_442_v2.0.window100kb.fasta

#bedtools getfasta -fi Meso_assembly/PvulgarisLaborOvalle_Meso_670_v1.0.fa -bed AM_window100kb.bed -fo PvulgarisLaborOvalle_Meso_670_v1.0.window100kb.fasta
