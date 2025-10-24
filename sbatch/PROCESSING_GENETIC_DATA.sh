# hg19 ref genome location:
## wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

#######################################
# GTEx - convert vcf to bfile
export TMPDIR=/scratch/jll1/tmp
cd /scratch/jll1/imQTL/genetic_data
cp /gpfs/data/pierce-lab/GTEx/GTEx_data/Genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz /scratch/jll1/imQTL/genetic_data

# annotate SNPs
input_vcf_file=/scratch/jll1/imQTL/genetic_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz
bcftools annotate ${input_vcf_file} --set-id '%CHROM:%POS:%REF:%ALT' | bcftools sort -Oz -o processed_gtex.vcf.gz
bcftools index -t processed_gtex.vcf.gz

# convert vcf to bfile
module load plink/2.0 
plink2 --vcf /scratch/jll1/imQTL/genetic_data/processed_gtex.vcf.gz --chr 1-22 --make-bed --out /scratch/jll1/imQTL/genetic_data/processed_gtex

###################
# below command writes out to TWAS project
###################
# module load plink/2.0
# plink2 --vcf processed_gtex.vcf.gz --chr 1-22 --make-bed --out /gpfs/data/huo-lab/BCAC/james.li/GTEx/genetic_data/processed_gtex

#######################################
# Identifying HEALS participants
R
library(data.table)
library(dplyr)
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_cov.RData")
IID_LIST <- gsub("X","",combined_cov$Sample_Name)
fam<-fread("/gpfs/data/phs/groups/Projects/GWAS/Bangladesh/genotype_data_imputed_multi_batches/Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr.fam") 
HEALS_FID_IID <- fam %>% filter(V2 %in% IID_LIST)
write.table(HEALS_FID_IID %>% select(V1,V2),file="/scratch/jll1/imQTL/genetic_data/HEALS_FID_IID.txt",quote=F,row.names=F,col.names=F,sep="\t")
R
q()
n

# filtering genetic data for HEALS participants
cd /scratch/jll1/imQTL/genetic_data
grep "D" /gpfs/data/phs/groups/Projects/GWAS/Bangladesh/genotype_data_imputed_multi_batches/Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005.bim | cut -f2 > remove_variants_Bangladesh.txt

# filtering participants and only retaining SNPs
module load plink/2.0
plink2 -bfile /gpfs/data/phs/groups/Projects/GWAS/Bangladesh/genotype_data_imputed_multi_batches/Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr \
--exclude /scratch/jll1/imQTL/genetic_data/remove_variants_Bangladesh.txt \
--keep /scratch/jll1/imQTL/genetic_data/HEALS_FID_IID.txt --recode vcf --output-chr MT --out /scratch/jll1/imQTL/genetic_data/Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr

# renaming chromosomes
rename_chr_file=/gpfs/data/pierce-lab/james.li/bcftools_chr_name_conv/chr_name_conv.txt
input_vcf_file=/scratch/jll1/imQTL/genetic_data/Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr.vcf
bcftools annotate --rename-chrs ${rename_chr_file} ${input_vcf_file} -Oz -o Bangladesh_rename_chr.vcf.gz

# running CrossMap
cd /scratch/jll1/imQTL/genetic_data/
fasta_file=/gpfs/data/pierce-lab/james.li/hg38/hg38.fa
input_chain=/gpfs/data/pierce-lab/james.li/liftOver/hg19ToHg38.over.chain.gz
CrossMap vcf ${input_chain} \
  /scratch/jll1/imQTL/genetic_data/Bangladesh_rename_chr.vcf.gz \
  ${fasta_file} \
  /scratch/jll1/imQTL/genetic_data/CrossMap_Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr.vcf

# converting to plink bfile
plink2 --vcf /scratch/jll1/imQTL/genetic_data/CrossMap_Bangladesh_HRC_hg37_imputed_r2gtp3_indels_mafgtp005_woAFdiff_wosexchr.vcf --allow-extra-chr --chr 1-22 --snps-only just-acgt --sort-vars --make-pgen --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw
plink2 -pfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw --chr 1-22 --make-bed --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw

######################################################
##### Obtaining final processed genotyping files #####
######################################################
# Examining overlap between Bangladesh and GTEx and renaming variants to make variant names consistent
R
library(dplyr)
library(data.table)
library(tidyr)

#####################################################
# importing gtex bim
gtex_bim <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/processed_gtex.bim")
gtex_bim <- gtex_bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":")) %>% mutate(newID=paste0("chr",ID))
# importing bangladesh bim
bangladesh_bim <- fread("/scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw.bim")
# checking number of intersecting records 
bangladesh_bim1 <- bangladesh_bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":"))
bangladesh_bim2 <- bangladesh_bim %>% mutate(ID=paste(V1,V4,V5,V6,sep=":"))
nrow(inner_join(gtex_bim,bangladesh_bim1,by=c("ID")))
nrow(inner_join(gtex_bim,bangladesh_bim2,by=c("ID")))

################# 
##### HEALS ##### 
################# 

#####################################################
# creating variant ID conversion dictionaries based A1 and A2 alleles designated by plink2
conversion_bangladesh <- bangladesh_bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":")) %>% mutate(newID=paste0("chr",ID)) %>% select(V2,newID)
gtex_bim <- fread("/scratch/jll1/imQTL/genetic_data/processed_gtex.bim")
conversion_gtex <- gtex_bim %>% mutate(ID=paste(V1,V4,V6,V5,sep=":")) %>% mutate(newID=paste0("chr",ID)) %>% select(V2,newID)

# writing out only Bangladesh dictionary since GTEx IDs are all identical between old and new IDs 
write.table(conversion_bangladesh,file="/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/HEALS_ID_CONVERT_DICT_A1A2.txt",quote=F,row.names=F,col.names=F,sep="\t")

# updating variant names for HEALS genotyping data
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw --update-name /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/HEALS_ID_CONVERT_DICT_A1A2.txt --make-bed --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename")
# filtering this file for genotype missingness and HWE
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename --geno 0.05 --hwe 1e-6 --make-bed --out /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename_geno0.05_hwe1e6")
# further filtering to exclude 4 duplicated variants (total 8 variants removed)
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/bfile/CrossMap_Bangladesh_raw_updatename_geno0.05_hwe1e6 --rm-dup exclude-mismatch --make-bed --out /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data")
# adding chr prefix, filtering MAF 0.01, selecting variants with no missingness
system("module load plink/2.0; plink2 -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data --maf 0.01 --geno 0 --output-chr chrM --keep-allele-order --make-pgen --out /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix")
# making an additional bfile copy 
system("module load plink/2.0; plink2 -pfile /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix --make-bed --out /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix_bfile")

################ 
##### GTEx ##### 
################ 
# filtering GTEx genetic data for genotype missingness and HWE
system("module load plink/2.0; plink2 -bfile /scratch/jll1/imQTL/genetic_data/processed_gtex --geno 0.05 --hwe 1e-6 --make-bed --out /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data")
# adding chr prefix and filtering MAF 0.01, selecting variants with no missingness
system("module load plink/2.0; plink2 -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data --maf 0.01 --geno 0 --output-chr chrM --keep-allele-order --make-pgen --out /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix")
# making an additional bfile copy 
system("module load plink/2.0; plink2 -pfile /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix --make-bed --out /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix_bfile")
q()
n
