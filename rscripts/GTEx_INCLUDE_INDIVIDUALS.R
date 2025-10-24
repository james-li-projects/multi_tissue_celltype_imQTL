library(data.table)
library(tidyverse)

# setting working directory and the directory to import DNAm samples from
setwd("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx")
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/methylation"

# initializing the loading function
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# importing genotyped samples
genotype.samplist <- (fread("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data.fam",header=F,sep="\t"))$V2

# initializing tissue list
tissue_list <- c("kidney","breast","colon","lung","prostate","wb","ovary")

# obtaining sample lists for each tissue type 
for (current_tissue in tissue_list) {
  print(paste("IDENTIFYING SAMPLES FOR TISSUE:",current_tissue))
  # importing DNAm samples
  noob_final_BMIQ <- loadRData(paste0(input_dir,"/noob_final_BMIQ_",current_tissue,"_2-6-2021.RData")) 

  # remove sample suffix so that we extract subject IDs
  colnames(noob_final_BMIQ)=str_extract(colnames(noob_final_BMIQ),'GTEX-\\w+')
  
  # load eligible GTEx sample lists
  load(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/eligible_gtex_samples/eligible_sample_list_",current_tissue,".RData"))
  
  # filtering DNAm data to only extract eligible samples for mQTL mapping
  dim(noob_final_BMIQ)
  filter_mat<-noob_final_BMIQ[,eligible_sample_list]
  dim(filter_mat)
  noob_final_BMIQ <- filter_mat
  dnam.samplist <- colnames(noob_final_BMIQ)
  
  # importing covariates, only import females for breast tissue
  if (current_tissue == "breast") {
    cov.samplist <- ((fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/covariates/",current_tissue,"_covariates_2-6-2021.csv"),header=T,sep=",")) %>% filter(SEX==2))$CollaboratorParticipantID
  } else {
    cov.samplist <- (fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/covariates/",current_tissue,"_covariates_2-6-2021.csv"),header=T,sep=","))$CollaboratorParticipantID
  }
  
  # identifying samples to use in mQTL analysis that have all the above data
  mQTL_sample_list <- intersect(intersect(genotype.samplist,dnam.samplist),cov.samplist)
  print(length(mQTL_sample_list))
  save(mQTL_sample_list, file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/mQTL_sample_list_",current_tissue,".RData"))
}
