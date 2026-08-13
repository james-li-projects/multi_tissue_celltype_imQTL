library(RNOmni)
library(stringr)
library(dplyr)
library(PCAForQTL)
library(Hmisc)
library(data.table)

input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/methylation"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS"
celltype_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac"

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
tissue_list <- c("wb")

# importing illumina annotations and making the backbone structure of a bed file
manifest<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/infinium-methylationepic-v-1-0-b5-manifest/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
manifest <- manifest %>% mutate(chr_string_hg19 = paste0("chr",CHR)) %>% filter(chr_string_hg19==CHR_hg38) %>% select(-chr_string_hg19)
bed_label <- manifest %>% select(CHR_hg38,Start_hg38,End_hg38,Name)

# make sure the chromosome value is prefixed with chr
colnames(bed_label) <- c("#chr","start","end","phenotype_id")

# performing the transform on every column
for (i in 1) {
  # importing DNAm data
  current_tissue <- tissue_list[i]
  print(paste("PERFORMING ANALYSIS FOR:",current_tissue))
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/mQTL_sample_list_",current_tissue,".RData"))

  # transforming beta values to m-values  
  print("Transforming beta values to M-values")
  beta_mat <- loadRData(paste0(input_dir,"/combined_beta_modid.RData"))
  beta_mat <- beta_mat[,colnames(beta_mat)%in%mQTL_sample_list]
  mval <- log2(beta_mat/(1 - beta_mat))

  # performing the inverse-normal transform on every row
  print(paste("Rank-based INT transforming"))
  transform_INT_beta_mat <- t(apply(X=mval,MARGIN=1, FUN=RankNorm))
  
  # generating and outputting DNAm bed file
  print(paste("Outputting DNAm .BED file"))
  transform_INT_beta_mat_df <- data.frame(transform_INT_beta_mat)
  colnames(transform_INT_beta_mat_df) <- gsub("X","",colnames(transform_INT_beta_mat_df))
  # arranging matrix by sample list order
  transform_INT_beta_mat_df <- transform_INT_beta_mat_df[,mQTL_sample_list]
  transform_INT_beta_mat_df$phenotype_id<-rownames(transform_INT_beta_mat_df)
  current_tissue_BED <- inner_join(bed_label,transform_INT_beta_mat_df, by = c("phenotype_id")) 
  # sorting bed file and removing cpgs on sex chromosomes
  current_tissue_BED <- current_tissue_BED %>% filter(!(`#chr`%in% c("chrX","chrY"))) %>% mutate(numeric_chr = as.numeric(gsub("chr","",`#chr`))) 
  current_tissue_BED <- current_tissue_BED %>% arrange(numeric_chr,start) %>% select(-numeric_chr)
  num_bed_samp <- ncol(current_tissue_BED)-4
  print(paste("Number samples in BED file:",num_bed_samp))
  # writing out bed file
  write.table(current_tissue_BED,file=paste0(output_dir,"/",current_tissue,".bed"),quote=F,row.names=F,col.names=T,sep="\t")
  
  # computing PCA for the current tissue DNAm matrix
  print(paste("Performing tissue-specific PCA"))
  prcompResult<-prcomp(t(transform_INT_beta_mat),center=TRUE,scale.=TRUE)
  PCs<-prcompResult$x
  # determining the appropriate number of PCs
  resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
  print(resultRunElbow)
  # extracting the appropriate number of PCs for the current tissue type
  QTL_PCS <- PCs[,1:resultRunElbow]
  save(QTL_PCS,file=paste0(output_dir,"/QTL_PCS_",current_tissue,".RData"))
  
  # loading cell type proportions for the current tissue
  load(paste0(celltype_dir,"/estF.",current_tissue,".RData"))
  estF.o <- data.frame(estF.o)
  estF.o$ID <- rownames(estF.o)
  QTL_PCS <- data.frame(QTL_PCS)
  QTL_PCS$ID <- rownames(QTL_PCS)
  
  # joining PCs with the cell type proportions
  joined_df <- inner_join(estF.o,QTL_PCS,by=c("ID")) %>% select(-ID)
  m <- as.matrix(joined_df)
  end_index <- head(which(grepl("PC",colnames(joined_df))),1)-1
  cormat<- rcorr(as.matrix(joined_df))$P[c(1:end_index),-c(1:end_index)]
  
  # outputting correlations between cell types and PCs
  save(cormat,file=paste0(celltype_dir,"/cormat_",current_tissue,".RData"))
}