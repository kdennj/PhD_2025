sessionInfo()
options(scipen=999)
library("GAPIT3")
dir <- paste0("/ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/bwa_Andean/gatk_filter_GWAS/GAPIT_chrs/pca3_bytrait")
dir
setwd(dir)

pheno <- read.csv(paste0("/ei/.project-scratch/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/KDJ_beans/bwa_Andean_M/gatk_filter_GWAS/GAPIT_chrs/23.06-phenotypes-GAPIT_all_in_numbers_Kate.txt"), header = T, sep = '\t') #header always true!

traits <- c("photo_sens", "Major_Seed_colour", "Seed_size", "Type", "Seed_coat", "Seed_pattern", "Sum_22", "Win23")

for(trait in traits){
print(trait)
phenocol<-grepl(trait, names(pheno)) 
mypheno<-cbind(pheno$TAXA, pheno[,phenocol])
print(names(mypheno))

dir <- paste0("/ei/projects/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/scratch/KDJ_beans/bwa_Andean/gatk_filter_GWAS/GAPIT_chrs/pca3_bytrait/", trait)
dir
dir.create(dir)
setwd(dir)

myGAPIT<-GAPIT(
file.G="Chr",
file.from=1,
file.to=11,
file.Ext.G="Andean_Haplotypecaller_filtered.vcf.biallelic_MAF01_thin5bp.vcf.gz.hapmap.txt_hmp.txt",
file.path="/ei/.project-scratch/3/3c493fe7-33a3-4d92-b57f-50d4e29c3391/KDJ_beans/bwa_Andean/gatk_filter_GWAS/GAPIT_chrs/",
Y=mypheno,
PCA.total = 3, 
Random.model=FALSE,
Major.allele.zero=T,
model = c("Blink", "FarmCPU", "MLM", "MLMM")
)
}