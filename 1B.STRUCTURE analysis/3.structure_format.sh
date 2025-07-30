#!/bin/bash -e
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -J jre
#SBATCH -p ei-short
#SBATCH --mem 20G
#SBATCH -c 1
#SBATCH -N 1

cd /KDJ_beans/bwa_Andean_M/gatk_filter_QUAL

#file= Andean_Haplotypecaller_filtered_maxhetero0.2_ID_LD-0.5_thin-w10.vcf
#outfile= STRUCTURE/Andean_MH0.2_LD0.5_w10.structure.txt

file=Andean_Haplotypecaller_filtered_maxhetero0.2_ID_LD-0.5_thin-w30.vcf
outfile=STRUCTURE-w30/Andean_MH0.2_LD0.5_w30.structure.txt

# convert to structure format
source jre-7.11

java -Xmx16g -Xms16g -jar /tgac/software/testing/pgdspider/2.1.1.2/x86_64/bin/PGDSpider2-cli.jar \
-inputfile ${file} \
-outputfile ${outfile} -spid STRUCTURE/vcfTOstructure.spid
