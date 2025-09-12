#!/bin/bash
#SBATCH -J popLDdecay --mem 200G
#SBATCH -N 1 # number of nodes
#SBATCH -t 00-04:00 # time (D-HH:MM)
#SBATCH -c 1 # number of tasks
#SBATCH -p ei-short
#SBATCH -o /KDJ_beans/bwa_Andean_M/gatk_filter_GWAS/popLDdecay/slurm.alleles.%N.%j.out # STDOUT
#SBATCH -e /KDJ_beans/bwa_Andean_M/gatk_filter_GWAS/popLDdecay/slurm.alleles.%N.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=kate.denning-james@earlham.ac.uk #send-to address

INPUT_VCF='/bwa_Andean_M/gatk_filter_GWAS/All_Andean_Haplotypecaller_GWAS_biallelic_MAF01_thin5bp.vcf'
#INPUT_VCF='/bwa_Andean_M/gatk_filter_QUAL/Andean_Haplotypecaller_filtered_maxhetero0.2_ID_LD-0.5_thin-w10.vcf'
#INPUT_VCF='/bwa_Andean_M/gatk_filter_GWAS/Chr01_Andean_Haplotypecaller_filtered.vcf.biallelic_MAF01_thin5bp.vcf.gz'


OUTPUT_PREFIX='bwa_Andean_M/gatk_filter_GWAS/popLDdecay/LDdecay_M'

Sub_population='bwa_Andean_M/gatk_filter_GWAS/popLDdecay/K2_Cluster_1_sample_uniq_names.csv'

#singularity exec ~/singularities/images/popLDdecay.img PopLDdecay -InVCF ${INPUT_VCF} -OutStat ${OUTPUT_PREFIX}

singularity exec ~/singularities/images/popLDdecay.img PopLDdecay -InVCF ${INPUT_VCF} -OutStat ${OUTPUT_PREFIX} -SubPop ${Sub_population}
