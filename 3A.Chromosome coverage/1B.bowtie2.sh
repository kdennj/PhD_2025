#!/bin/bash
#SBATCH --mail-type=END, FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address
#SBATCH -p ei-long
#SBATCH -c 4
#SBATCH --mem 100G
#SBATCH -J bowtie2

source bowtie-2.4.5
source perl-5.22.1
source samtools-1.7

cd /KDJ_beans/

DIR=/KDJ_beans/trim_reads
OUTDIR=/KDJ_beans/bowtie2_meso_andean

#mapfile -t sample < /KDJ_beans/samples_withoutpractise.txt

sample=$1

bowtie2 -q -N 0 -X 450 --fr --no-unal --no-mixed --no-discordant --no-contain --no-overlap \
-x ref_genomes/Andean_Meso_ref_genome/Andean_Meso_ref_bowtie \
-1 <(zcat ${DIR}/${sample}_L001_R1_val_1.fq.gz ${DIR}/${sample}_L002_R1_val_1.fq.gz) \
-2 <(zcat ${DIR}/${sample}_L001_R2_val_2.fq.gz ${DIR}/${sample}_L002_R2_val_2.fq.gz) \
-S - | samtools view -bS -o ${OUTDIR}/${sample}_bowtie2.bam && touch ${sample}.bowtie2done 


#-X 450 means maximum insert size of 450 - to allow margin for indels 
#-fr means if a candidate paired-end alignment where mate 1 appears upstream of reverse compliment of mate2 and insert length constraints are met then alignment is valid
#--no-unal suppresses sam records for reads that failed to align 
#-N Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1
#--no-mixed By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior
#--no-discordant bowtie2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior.
#--no-contain If one mate alignment contains the other, consider that to be non-concordant
#--no-overlap If one mate alignment overlaps the other at all, consider that to be non-concordant

#The --strata and --best options do not apply in paired-end mode.
#-- best means reported singleton alignments are best in terms of stratum (number mismatches) 
#-- strata means If many valid arguments exist but aren't reportable and they fall into one or more alignment stratum then only report the best stratum 