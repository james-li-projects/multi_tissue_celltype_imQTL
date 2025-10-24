library(data.table)
library(tidyverse)
library(EpiSCORE)
library(reshape)
library(EpiDISH)

meth_input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/methylation"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac"

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
tissue_list <- c("kidney","breast","colon","lung","prostate","wb","ovary")

# loading centEpiFibEndoIC.m
load("/gpfs/data/pierce-lab/james.li/imQTL/teschendorff_reference/centEpiFibEndoIC.Rd")

tissue_mref_list <- vector(length=length(tissue_list),mode="list")
tissue_mref_list[[1]]<-EpiSCORE::Kidney_Mref.m
tissue_mref_list[[2]]<-EpiSCORE::mrefBreast.m
tissue_mref_list[[3]]<-EpiSCORE::Colon_Mref.m
tissue_mref_list[[4]]<-EpiSCORE::mrefLung.m
tissue_mref_list[[5]]<-EpiSCORE::mrefProstate.m
tissue_mref_list[[6]]<-EpiDISH::centDHSbloodDMC.m
tissue_mref_list[[7]]<-centEpiFibEndoIC.m

# estimating cell type fractions for all tissues using EpiSCORE except whole blood and ovary
for (i in 1:5) {
  # identifying current tissue type
  current_tissue <- tissue_list[i]
  current_tissue_mref <- tissue_mref_list[[i]]
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/mQTL_sample_list_",current_tissue,".RData"))
  
  print(paste("Inferring cell type fractions for",current_tissue))
  # restrict samples 
  noob_final_BMIQ <- loadRData(paste0(meth_input_dir,"/noob_final_BMIQ_",current_tissue,"_2-6-2021.RData"))
  colnames(noob_final_BMIQ)=str_extract(colnames(noob_final_BMIQ),'GTEX-\\w+')
  noob_final_BMIQ <- noob_final_BMIQ[,colnames(noob_final_BMIQ)%in%mQTL_sample_list]
  
  # inferring cell type fractions
  constAvBetaTSS.o <-constAvBetaTSS(noob_final_BMIQ,type='850k')
  wRPC.o <- wRPC(constAvBetaTSS.o,ref=current_tissue_mref,useW=TRUE,wth=0.4,maxit=200)
  estF.o <- wRPC.o$estF
  #melt data for plotting
  difficulty_data <- estF.o %>%
    t() %>%
    as.data.frame %>% 
    t() 
  plot_df <- as.data.frame.table(difficulty_data)
  colnames(plot_df) <- c("Sample_ID","Cell Type","Proportion")
  png(paste0(output_dir,"/cell_type_fraction_",current_tissue,".png"),units="in",height=4,width=4,res=1200)
  p <- ggplot(plot_df, aes(x=`Cell Type`, y=Proportion)) + geom_boxplot() + theme_classic()
  print(p)
  dev.off()
  # saving estimated cell type fractions
  save(estF.o,file=paste0(output_dir,"/estF.",current_tissue,".RData"))
}

# estimating cell type fractions for whole blood and ovary
for (i in 6:7) {
  current_tissue <- tissue_list[i]
  current_tissue_mref <- tissue_mref_list[[i]]
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/mQTL_sample_list_",current_tissue,".RData"))
  print(paste("Inferring cell type fractions for",current_tissue))
  noob_final_BMIQ <- loadRData(paste0(meth_input_dir,"/noob_final_BMIQ_",current_tissue,"_2-6-2021.RData"))
  colnames(noob_final_BMIQ)=str_extract(colnames(noob_final_BMIQ),'GTEX-\\w+')
  noob_final_BMIQ <- noob_final_BMIQ[,colnames(noob_final_BMIQ)%in%mQTL_sample_list]
  wRPC.o <- epidish(beta.m = noob_final_BMIQ, ref.m = current_tissue_mref, method = 'RPC')
  estF.o <- wRPC.o$estF
  #melt data for plotting
  difficulty_data <- estF.o %>%
    t() %>%
    as.data.frame %>% 
    t() 
  plot_df <- as.data.frame.table(difficulty_data)
  colnames(plot_df) <- c("Sample_ID","Cell Type","Proportion")
  png(paste0(output_dir,"/cell_type_fraction_",current_tissue,".png"),units="in",height=4,width=4,res=1200)
  p <- ggplot(plot_df, aes(x=`Cell Type`, y=Proportion)) + geom_boxplot() + theme_classic()
  print(p)
  dev.off()
  # saving estimated cell type fractinos
  save(estF.o,file=paste0(output_dir,"/estF.",current_tissue,".RData"))
}
