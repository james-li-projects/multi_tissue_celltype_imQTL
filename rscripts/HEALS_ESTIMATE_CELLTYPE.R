library(data.table)
library(tidyverse)
library(EpiSCORE)
library(reshape)
library(EpiDISH)

meth_input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/methylation"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac"

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
tissue_list <- c("wb")

tissue_mref_list <- vector(length=length(tissue_list),mode="list")
tissue_mref_list[[1]]<-EpiDISH::centDHSbloodDMC.m

# estimating cell type fractions for whole blood
for (i in 1) {
  current_tissue <- tissue_list[i]
  current_tissue_mref <- tissue_mref_list[[i]]
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/mQTL_sample_list_",current_tissue,".RData"))
  print(paste("Inferring cell type fractions for",current_tissue))
  beta_mat <- loadRData(paste0(meth_input_dir,"/combined_beta_modid.RData"))
  beta_mat <- beta_mat[,colnames(beta_mat)%in%mQTL_sample_list]
  wRPC.o <- epidish(beta.m = beta_mat, ref.m = current_tissue_mref, method = 'RPC')
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
