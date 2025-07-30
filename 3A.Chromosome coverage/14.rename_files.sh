#!/bin/bash -e


cd /KDJ_beans/bowtie2_meso_andean/unique/qualimap

for file in $(cat /KDJ_beans/samples_uniq_all.txt);
do
  cd ${file}_unique.bam
  mv genome_results.txt ${file}_genome_results.txt
  mv ${file}_genome_results.txt  ../all_qual_res
  cd /KDJ_beans/bowtie2_meso_andean/unique/qualimap

done
