#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -c 1 --mem 50G 
#SBATCH -J bamtools
#SBATCH --mail-user=denning@nbi.ac.uk
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --array=1-145
#SBATCH -t 1-00:00

source bamtools-2.5.1_CBG

mapfile -t file < samples_uniq_all.txt

bamtools index -in ${file[${SLURM_ARRAY_TASK_ID}]}_A.mdup.bam