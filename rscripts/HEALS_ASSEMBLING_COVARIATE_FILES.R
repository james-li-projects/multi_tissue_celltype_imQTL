library(data.table)
library(dplyr)
library(stringr)
library(RNOmni)

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS"

tissue_list <- c("wb")

for (current_tissue in tissue_list) {
  print(paste("PROCESSING COVARIATE FILE FOR TISSUE:",current_tissue))
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/mQTL_sample_list_",current_tissue,".RData"))
  
  # importing HEALS covariates
  current_tissue_cov <- loadRData(paste0(input_dir,"/covariates/combined_cov_modid.RData"))
  current_tissue_cov <- current_tissue_cov %>% select(Sample_Name,Sex) %>% filter(Sample_Name%in%mQTL_sample_list)
  print(paste("TOTAL NUMBER OF RETAINED SAMPLES FOR mQTL ANALYSIS:",nrow(current_tissue_cov)))
  
  # importing genotyping PCs
  genotyping_PC <- fread(paste0(input_dir,"/genetic_data/PLINK2_PCA_30.eigenvec"),header=T) %>% select(IID, PC1, PC2, PC3, PC4, PC5)
  colnames(genotyping_PC)[1:6] <- c("Sample_Name","G_PC1","G_PC2","G_PC3","G_PC4","G_PC5")
  genotyping_PC <- genotyping_PC %>% select(Sample_Name,G_PC1,G_PC2,G_PC3,G_PC4,G_PC5)
  # importing tissue-specific PCs
  current_tissue_PC <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/QTL_PCS_",current_tissue,".RData"))
  current_tissue_PC_df <- data.frame(current_tissue_PC)
  current_tissue_PC_df$Sample_Name <- rownames(current_tissue_PC_df)
  # assembling and outputing covariate files
  combined_cov <- inner_join(inner_join(current_tissue_cov,genotyping_PC,by=c("Sample_Name")),current_tissue_PC_df,by=c("Sample_Name"))
  # if the tissue is prostate or ovary, exclude the sex variable since all samples are from males
  # # # print(table(current_tissue_cov$Sex))
  if (current_tissue %in% c("prostate","ovary","breast")) {
    combined_cov <- combined_cov %>% select(-Sex)
  } else {
    print("Tissue does not only have samples that are sex-specific (e.g. prostate, ovary, or breast)")
  }
  
  # keeping the order of the intersecting sample list
  combined_cov <- data.frame(combined_cov)
  rownames(combined_cov) <- combined_cov$Sample_Name
  combined_cov <- combined_cov[mQTL_sample_list,]
  
  # further manipulating the covariate data.frame to match outputs for tensorQTL
  t_combined_cov <- t(combined_cov)
  colnames(t_combined_cov) <- t_combined_cov[1, ]
  t_combined_cov <- t_combined_cov[-1, ]
  write.table(t_combined_cov,file=paste0(output_dir,"/processed_covariates_",current_tissue,".txt"),quote=F,row.names=T,col.names=T,sep="\t")
  
  print(paste("PROCESSING CELL TYPE FILE FOR TISSUE:",current_tissue))
  load(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac/estF.",current_tissue,".RData"))
  estF.o <- data.frame(estF.o)
  estF.o <- estF.o[mQTL_sample_list,]
  # performing the ranked-based inverse normal transform (INT) on every cell type fraction column
  estF.o <- (apply(X=estF.o,MARGIN=2, FUN=RankNorm))
  
  # writing out cell type fractions for interaction analysis 
  for (CT_index in seq(colnames(estF.o))) {
    current_CT <- colnames(estF.o)[CT_index]
    current_tissue_CT_interaction_df <- data.frame(rownames(estF.o),estF.o[,CT_index])
    write.table(current_tissue_CT_interaction_df,file=paste0(output_dir,"/processed_interactions_",current_tissue,"_",current_CT,".txt"),quote=F,row.names=F,col.names=F,sep="\t")
  }
}