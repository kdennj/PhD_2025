#!/bin/bash
#SBATCH -p ei-long # partition (queue)
#SBATCH -c 1
#SBATCH -t 1
#SBATCH --mem 100G # memory pool for all cores
#SBATCH -t 4-00:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk # send-to address
#SBATCH --array=0-1

#1> tells the computer to kindly pipe the standard output into whatever file you tell it.
source Structure-2.3.5

K=(5 9)

cd /KDJ_beans/bwa_Andean_M/gatk_filter_QUAL/STRUCTURE-w30

for ((i=61; i<=80; i++))
do
		echo Beginning loop K=${K[${SLURM_ARRAY_TASK_ID}]}_${i}
		structure -K ${K[${SLURM_ARRAY_TASK_ID}]} -o Results/K${K[${SLURM_ARRAY_TASK_ID}]}_R${i} 1> ${K[${SLURM_ARRAY_TASK_ID}]}_${i}_structure.log
done
