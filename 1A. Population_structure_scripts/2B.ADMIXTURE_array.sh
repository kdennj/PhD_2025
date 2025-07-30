#!/bin/bash
#SBATCH -p ei-medium # partition (queue)
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem 30G # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk # send-to address
#SBATCH --array=0-9

source admixture-1.3.0

K=(2 3 4 5 6 7 8 9 10)

cd ADMIXTURE-cluster_1/Results

BED_FILE=Andean_Haplotypecaller_filtered_maxhetero0.2_ID_LD-0.5_C1_MAF10_thin-w10
OUTDIR=ADMIXTURE

for ((i=1; i<=10; i++))
do
		echo Beginning loop K=${K[${SLURM_ARRAY_TASK_ID}]}_${i}
    mkdir K${K[${SLURM_ARRAY_TASK_ID}]}_R${i}_$OUTDIR
    cd K${K[${SLURM_ARRAY_TASK_ID}]}_R${i}_$OUTDIR
    admixture --cv ../$BED_FILE.bed ${K[${SLURM_ARRAY_TASK_ID}]}
    cd ..
done
