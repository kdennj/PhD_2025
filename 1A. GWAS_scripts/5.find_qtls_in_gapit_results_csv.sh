#!/bin/bash
#SBATCH -p ei-short -c 1 --mem 1G

source python-3.5.1 && source bedtools-2.28.0

windows=/hpc-home/denning/KDJ_Scratch_project/KDJ_beans/ref_genomes/Andean_Meso_ref_genome/Andean_genome_10kb_sliding_winds.txt

input=$1

zcat ${input} | grep -v 'Pos' | tr ',' '\t' | awk -F '\t' '{print $2 "\t" $3 "\t" $3+1 "\t" log($4)/log(10)*-1}' > ${input}.bed4
bedtools map -a ${windows} -b <(sort -k1,1 -k2,2n ${input}.bed4) -c 4,4 -o max,collapse | awk '{if ($4 >= 3) print $1 "_" $2 "\t" $5}' > ${input}.bed4.maxover5.txt

cat ${input}.bed4.maxover5.txt | awk -F '\t' -v input_file=${input} '{
text_str = $1;
split($2, num_values, ",");
num_ge_9 = num_ge_8 = num_ge_7 = num_ge_6 = num_ge_5 = num_ge_4 = num_ge_3 = 0;
for (i in num_values) {
if (num_values[i] >= 9) num_ge_9++;
if (num_values[i] >= 8) num_ge_8++;
if (num_values[i] >= 7) num_ge_7++;
if (num_values[i] >= 6) num_ge_6++;
if (num_values[i] >= 5) num_ge_5++;
if (num_values[i] >= 4) num_ge_4++;
if (num_values[i] >= 3) num_ge_3++;
}
print input_file "\t" text_str "\t" num_ge_9 "\t" num_ge_8 "\t" num_ge_7 "\t" num_ge_6 "\t" num_ge_5 "\t" num_ge_4 "\t" num_ge_3;
}' > ${input}.bed4.maxover5_output.txt