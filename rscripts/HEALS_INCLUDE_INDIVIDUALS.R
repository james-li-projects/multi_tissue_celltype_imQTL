library(data.table)
library(tidyverse)
library(tidyr)

# setting directory to import DNAm samples from
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS"

# initializing the loading function
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# importing genotyped samples
genotype_fam <- (fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data.fam",header=F,sep="\t")) %>% separate(V2, into=c("FID","IID"),remove=F) %>% mutate(Sample_Name=paste0("X",IID))
genotype.samplist <- genotype_fam$Sample_Name

# obtaining sample list 
current_tissue="wb"
print(paste("IDENTIFYING SAMPLES FOR TISSUE:",current_tissue))
# importing DNAm samples and extracting sample names
beta_val <- loadRData(paste0(input_dir,"/methylation/combined_beta.RData")) 
dnam.samplist <- colnames(beta_val)
# importing covariates
imported_cov <- loadRData(paste0(input_dir,"/covariates/combined_cov.RData"))
cov.samplist <- imported_cov$Sample_Name
# identifying samples to use in mQTL analysis that have all the above data
mQTL_sample_list_XID <- intersect(intersect(genotype.samplist,dnam.samplist),cov.samplist)
print(length(mQTL_sample_list_XID))

# ID conversion dictionary
dict <- genotype_fam %>% select(V2,Sample_Name) %>% filter(Sample_Name %in% mQTL_sample_list_XID)

# filtering covariate file for intersecting samples and giving them IDs from the genotype file
rownames(imported_cov) <- imported_cov$Sample_Name
imported_cov <- imported_cov[dict$Sample_Name,]
combined_cov_modid <- inner_join(imported_cov,dict,by=c("Sample_Name"))
combined_cov_modid <- combined_cov_modid %>% mutate(Sample_Name=V2) %>% select(-V2)
save(combined_cov_modid, file = paste0(input_dir,"/covariates/combined_cov_modid.RData"))

# filtering methylation matrix for intersecting samples and giving them IDs from the genotype file
combined_beta_modid<-beta_val[,dict$Sample_Name]
colnames(combined_beta_modid) <- dict$V2
save(combined_beta_modid, file = paste0(input_dir,"/methylation/combined_beta_modid.RData"))

# saving RData object containing intersecting IDs
mQTL_sample_list <- dict$V2
save(mQTL_sample_list, file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/mQTL_sample_list_",current_tissue,".RData"))
