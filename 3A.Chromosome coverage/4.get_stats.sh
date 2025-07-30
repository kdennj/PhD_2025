for i in *bam;

do

sbatch -J stats --mem 20G -c 1 -N 1 -p ei-medium -t 1-0 --wrap "source samtools-1.7;samtools flagstat ${i} > stats.${i}.txt";

done
