#!/bin/bash -e

##FOR LOOP 2 TO 10 (repeat 10 times):
source admixture-1.3.0

#admixture –cv output.bed ${i}
BED_FILE=Andean_Haplotypecaller_filtered_maxhetero0.2_ID_LD-0.5_thin-w10.bed.bed
OUTDIR=ADMIXTURE
for ((K=2; K<=10; K+=1)); do
echo Starting run for K=${K}
for ((i=1; i<=10; i++)); do
echo Beginning loop K="$K"_"$i"
mkdir K${K}_R${i}_$OUTDIR
cd K${K}_R${i}_$OUTDIR
sbatch -J admx -o admix$K$i.log --mem 30G -c 1 -N 1 -p ei-medium --wrap "source admixture-1.3.0;admixture --cv ../$BED_FILE ${K}"
cd ..
done
done
