#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -J minimap
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH --mem 20G # memory pool for all cores
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address


#source minimap2-2.7
#source minimap2-2.22_CBG

source minimap2-2.24

cd /hpc-home/denning/KDJ_Scratch_project/KDJ_beans/ref_genomes


minimap2 -x asm10 Meso_assembly/PvulgarisLaborOvalle_Meso_670_v1.0.fa Pvulgaris_Andean_442_v2.0.window100kb.fasta > Andean_Meso_ref_genome/A_100kbwindows-over-M.paf

#minimap2 -x asm10 Andean_assembly/Pvulgaris_Andean_442_v2.0.fa PvulgarisLaborOvalle_Meso_670_v1.0.window100kb.fasta > M_100kbwindows-over-A.paf



###Notes
#-x asm10:
#comparing two genome assemblies that are about 90% or more identical — perfect for intra-species comparisons
#Andean_assembly/Pvulgaris_Andean_442_v2.0.fa:
#The reference genome you're aligning to 
#PvulgarisLaborOvalle_Meso_670_v1.0.window100kb.fasta:
#These are the query sequences
#-o M_100kbwindows-over-A.paf:
#.paf file (Pairwise mApping Format), which stores the alignments.