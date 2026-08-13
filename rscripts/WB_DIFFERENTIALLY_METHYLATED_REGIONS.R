# initializing packages
library(data.table)
library(tidyr)
library(dplyr)
library(tidyverse)
library(sva)
library(mgcv)
library(nlme)
library(genefilter)
library(BiocParallel)
library(ggplot2)
library(RNOmni)
library(stringr)
library(PCAForQTL)
library(Hmisc)
library(gap)

#####################
# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#####################

#########
# HEALS #
#########
# importing methylation and covariate data
#load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_cov.RData")
load("/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/data/combined_mval.RData")
noob_final_BMIQ <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/methylation/noob_final_BMIQ_wb_2-6-2021.RData"))

heals_cpgs <- rownames(combined_mval)
gtex_cpgs <- rownames(noob_final_BMIQ)
pan_dataset_cpgs <- intersect(heals_cpgs,gtex_cpgs)


############################################
# importing methylation and covariate data #
############################################
# loading sample list for gtex blood
mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/mQTL_sample_list_wb.RData"))

# transforming beta values to m-values  
print("Transforming beta values to M-values")
noob_final_BMIQ <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/methylation/noob_final_BMIQ_wb_2-6-2021.RData"))
colnames(noob_final_BMIQ)=str_extract(colnames(noob_final_BMIQ),'GTEX-\\w+')
noob_final_BMIQ <- noob_final_BMIQ[,colnames(noob_final_BMIQ)%in%mQTL_sample_list]
mval <- log2(noob_final_BMIQ/(1 - noob_final_BMIQ))
mval <- mval[pan_dataset_cpgs,]

# importing GTEx covariates
cov <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/covariates/wb_covariates_2-6-2021.csv"),header=T,sep=",")
cov <- cov %>% select(CollaboratorParticipantID,SEX,AGE,Sample_Plate) %>% filter(CollaboratorParticipantID%in%mQTL_sample_list) %>% mutate(SEX = SEX-1) %>% rename(Sample_Name=CollaboratorParticipantID)
print(paste("TOTAL NUMBER OF RETAINED SAMPLES FOR DiffMeth ANALYSIS:",nrow(cov)))

# initializing table to store all results
all_sig_results <- data.frame()

# importing cell type scores
for (current_CT_name in c("B","CD8T","CD4T","Neutro","Mono","NK")) {
  print(current_CT_name)
  CT <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/processed_interactions_wb_",current_CT_name,".txt")) %>% mutate(Sample_Name=V1,CT_SCORE=V2)
  
  # reading resulting covariates df
  tmp_cov <- inner_join(cov,CT,by=c("Sample_Name"))
  
  # making sure mval and cov files are aligned
  mval <- mval[,tmp_cov$Sample_Name]
  
  # EWAS MVAL 
  fmla <- as.formula(paste0("~ tmp_cov$CT_SCORE + factor(tmp_cov$SEX) + tmp_cov$AGE + factor(tmp_cov$Sample_Plate)"))
  design <- model.matrix(fmla)
  library(limma)
  fit <- lmFit(mval, design)
  fit <- eBayes(fit)
  results <- topTable(fit, n=dim(mval)[1], sort.by="P",adjust="BH", coef=paste0("tmp_cov$CT_SCORE"), confint=TRUE)

  # plotting QQ-plot
  observed_pvalues <- na.omit(results$P.Value)
  file_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/DiffMeth/qqplot_", current_CT_name, ".png")
  
  lambda <- round(median(qchisq(1 - observed_pvalues, df = 1)) / qchisq(0.5, df = 1), 3)
  
  png(file = file_path)
  qqunif(observed_pvalues, main = "QQ Plot of P-values")
  text(x = 0.3, y = max(-log10(observed_pvalues)) * 0.9, labels = bquote(lambda == .(lambda)), col = "red", cex = 1.2)
  dev.off()
  
  # filtering for significant results
  sig_results<-results%>%filter(adj.P.Val<0.05) %>% mutate(Name=rownames(.),CellType=current_CT_name)
  
  all_sig_results <- rbind(all_sig_results,sig_results)
}

#######################
# importing wb imQTLs #
#######################
setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc"))
file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
file_list <- file_list[!grepl("Eosino",file_list)]
imQTL_cpg_df <- data.frame()
for (combination in file_list) {
  print(combination)
  tmp_df <- fread(combination)
  
  # filtering for bonferroni significant imQTLs
  tmp_df <- tmp_df %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni"))
  tmp_df <- tmp_df %>% filter(pval_adj_bonf < 0.05)

  # assembling a DF containing all the cpgs
  tmp_imQTL_cpg_df <- tmp_df %>% select(phenotype_id,b_g,b_gi) %>% mutate(combination=gsub("tensorQTL_imQTL_wb_","",gsub(".cis_qtl_top_assoc.txt.gz","",combination)))
  imQTL_cpg_df <- rbind(imQTL_cpg_df,tmp_imQTL_cpg_df)
}

######################
# importing wb mQTLs #
######################
setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/mQTL/top_assoc"))
mQTL_cpg_df <- fread("tensorQTL_mQTL_wb.cis_qtl.txt.gz") %>% mutate(pval_adj_bh=p.adjust(pval_beta,method="fdr")) %>% filter(pval_adj_bh < 0.05) %>% select(phenotype_id)

###############################
# identifying all differentially methylated CpGs
DiffMeth_CpG_list <- unique(all_sig_results$Name)

# examining proportion of imQTL and mQTL CpGs that are located in DMS
length(intersect(DiffMeth_CpG_list,unique(imQTL_cpg_df$phenotype_id)))
length(unique(imQTL_cpg_df$phenotype_id))

length(intersect(DiffMeth_CpG_list,unique(mQTL_cpg_df$phenotype_id)))
length(unique(mQTL_cpg_df$phenotype_id))

# computing the proportion of mQTLs in DMS
mQTL_DiffMeth_prop <- length(intersect(DiffMeth_CpG_list,unique(mQTL_cpg_df$phenotype_id))) / length(unique(mQTL_cpg_df$phenotype_id))


############################################
all_ct_df <- data.frame()

n_base_1<-length(DiffMeth_CpG_list)
n_base_2<-length(unique(pan_dataset_cpgs))

n_mqtl_1<-length(intersect(DiffMeth_CpG_list,unique(mQTL_cpg_df$phenotype_id)))
n_mqtl_2<-length(unique(mQTL_cpg_df$phenotype_id))

n_imqtl_1<-length(intersect(DiffMeth_CpG_list,unique(imQTL_cpg_df$phenotype_id)))
n_imqtl_2<-length(unique(imQTL_cpg_df$phenotype_id))

#n_a_1<-length(intersect(DiffMeth_CpG_list,unique((imQTL_cpg_df %>% filter(sign(b_g)==sign(b_gi)))$phenotype_id)))
#n_a_2<-length(unique((imQTL_cpg_df %>% filter(sign(b_g)==sign(b_gi)))$phenotype_id))

#n_b_1<-length(intersect(DiffMeth_CpG_list,unique((imQTL_cpg_df %>% filter(sign(b_g)!=sign(b_gi)))$phenotype_id)))
#n_b_2<-length(unique((imQTL_cpg_df %>% filter(sign(b_g)!=sign(b_gi)))$phenotype_id))

all_ct_df <- rbind(all_ct_df,c(n_base_1,n_base_2,"EPIC array"))
all_ct_df <- rbind(all_ct_df,c(n_mqtl_1,n_mqtl_2,"mQTL"))
all_ct_df <- rbind(all_ct_df,c(n_imqtl_1,n_imqtl_2,"imQTL"))
#all_ct_df <- rbind(all_ct_df,c(n_a_1,n_a_2,"Same Direction"))
#all_ct_df <- rbind(all_ct_df,c(n_b_1,n_b_2,"Different Direction"))

colnames(all_ct_df) <- c("DiffMeth","TotalCpG","Classification")
all_ct_df <- all_ct_df %>% mutate(DiffMeth=as.numeric(DiffMeth),TotalCpG=as.numeric(TotalCpG)) %>% mutate(Prop_DiffMeth=(DiffMeth)/(TotalCpG),Prop_DiffMeth_SE=sqrt(Prop_DiffMeth*(1-Prop_DiffMeth)/TotalCpG),Prop_DiffMeth_L=Prop_DiffMeth-1.96*Prop_DiffMeth_SE,Prop_DiffMeth_U=Prop_DiffMeth+1.96*Prop_DiffMeth_SE)

# Ensure Classification is a factor with the specified order
all_ct_df$Classification <- factor(all_ct_df$Classification, levels = c(
  "EPIC array",
  "mQTL",
  "imQTL"#,
  #"Same Direction",
  #"Different Direction"
))

library(dplyr)
library(ggplot2)

# Add a column with formatted percent labels
all_ct_df <- all_ct_df %>%
  mutate(percent_label = paste0(round(Prop_DiffMeth * 100, 1), "%"))

# write tsv of these values
write.table(all_ct_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/DiffMeth/counts.tsv")


imqtl_count_subset <- subset(all_ct_df, Classification == "imQTL")
mqtl_count_subset  <- subset(all_ct_df, Classification == "mQTL")

# Two-proportion test
prop_test <- prop.test(
  x = c(imqtl_count_subset$DiffMeth, mqtl_count_subset$DiffMeth),
  n = c(imqtl_count_subset$TotalCpG, mqtl_count_subset$TotalCpG),
  alternative = "two.sided"   # tests whether imQTL > mQTL
)

# Create the plot
p <- ggplot(all_ct_df, aes(x = Classification, y = Prop_DiffMeth, fill = Classification)) +
  geom_col(color = "black", width = 0.7) +  # Add bar chart
  geom_errorbar(aes(ymin = Prop_DiffMeth_L, ymax = Prop_DiffMeth_U), width = 0.2) +  # Overlay error bars
  geom_point(size = 2, color = "black") +  # Optional: highlight point estimate
  geom_text(aes(label = percent_label), vjust = -2, size = 7) +  # Add percent label above bar
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.ticks.x = element_line(size = 0.8),
    legend.position = "none"
  ) +
  labs(y = "% of CpGs in sites associated with\ninferred cell type proportions") +
  ylim(0,0.65)

# Save the plot
ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/DiffMeth/barplot_errorbar_allimqtl.png",
  plot = p,
  width = 6, height = 6, dpi = 1200
)

