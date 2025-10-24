library(data.table)
library(dplyr)
library(tidyr)

#######################################
# extracting numbers of mQTLs
mQTL_feature_df <- data.frame()
mQTL_cpg_df <- data.frame()

for (Dataset in c("GTEx","HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/mQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl.txt.gz",list.files())]
  for (file_name in file_list) {
    print(file_name)
    tmp_df <- fread(file_name)
    tmp_df <- tmp_df %>% mutate(pval_adj_bh=p.adjust(pval_beta,method="fdr")) %>% filter(pval_adj_bh < 0.05)
    num_mQTL <- nrow(tmp_df)
    print(num_mQTL)
    tmp_mQTL_feature_df <- data.frame(Dataset,file_name,num_mQTL) %>% mutate(Tissue=gsub(".cis_qtl.txt.gz","",file_name)) %>% mutate(Tissue=gsub("tensorQTL_mQTL_","",Tissue))
    # assembling feature DF 
    mQTL_feature_df <- rbind(mQTL_feature_df,tmp_mQTL_feature_df)
    
    # assembling a DF containing all the cpgs
    tmp_mQTL_cpg_df <- tmp_df %>% select(phenotype_id) %>% mutate(file_name=file_name) %>% mutate(Tissue=gsub(".cis_qtl.txt.gz","",file_name)) %>% mutate(Tissue=gsub("tensorQTL_mQTL_","",Tissue)) %>% mutate(Dataset=Dataset)
    mQTL_cpg_df <- rbind(mQTL_cpg_df,tmp_mQTL_cpg_df)
  }
}

#######################################
# making the tissue names in DFs more beautiful
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="breast"] <- "Breast"
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="colon"] <- "Colon"
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="lung"] <- "Lung"
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="kidney"] <- "Kidney"
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="prostate"] <- "Prostate"
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="wb"] <- "Whole Blood"
mQTL_feature_df$Tissue[mQTL_feature_df$Tissue=="ovary"] <- "Ovary"
# making the tissue names in the mQTL_feature_df more beautiful
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="breast"] <- "Breast"
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="colon"] <- "Colon"
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="lung"] <- "Lung"
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="kidney"] <- "Kidney"
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="prostate"] <- "Prostate"
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="wb"] <- "Whole Blood"
mQTL_cpg_df$Tissue[mQTL_cpg_df$Tissue=="ovary"] <- "Ovary"

#######################################
# making copies of the feature and cpg dfs of mQTLs before filtering out tissues with no imQTLs
copy_mQTL_feature_df <- mQTL_feature_df
copy_mQTL_cpg_df <- mQTL_cpg_df
# filtering out tissue/celltype combinations with no imQTLs
mQTL_feature_df<-mQTL_feature_df %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Whole Blood")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Breast")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Kidney"))
mQTL_cpg_df<-mQTL_cpg_df %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Whole Blood")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Breast")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Kidney"))

#######################################
# saving DFs
save(mQTL_cpg_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_cpg_df.RData")
save(mQTL_feature_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_feature_df.RData")


################################
# Define the mapping for the Tissue column
tissue_mapping <- c(
  "Lung" = "Lung (n=190)",
  "Colon" = "Colon (n=189)",
  "Ovary" = "Ovary (n=140)",
  "Prostate" = "Prostate (n=105)",
  "Kidney" = "Kidney (n=47)",
  "Breast" = "Breast (n=36)",
  "Whole Blood" = "Whole Blood (n=1,182)"
)

# Apply the mapping to the Tissue column
modified_mQTL_feature_df <- copy_mQTL_feature_df
modified_mQTL_feature_df$Tissue <- recode(modified_mQTL_feature_df$Tissue, !!!tissue_mapping)
modified_mQTL_feature_df <- modified_mQTL_feature_df %>% mutate(Tissue=ifelse(Tissue=="Whole Blood (n=1,182)"&Dataset=="GTEx","Whole Blood (n=47)",Tissue))
# Define custom colors for specific Tissue values
tissue_colors <- c(
  "Breast (n=36)" = "turquoise3",
  "Lung (n=190)" = "yellowgreen",
  "Colon (n=189)" = "sienna",
  "Ovary (n=140)" = "pink3",
  "Prostate (n=105)" = "lightgray",
  "Whole Blood (n=1,182)" = "magenta3",
  "Whole Blood (n=47)" = "magenta3",
  "Kidney (n=47)" = "turquoise2"
)

####################
# generating plots #
####################
library(ggplot2)
library(scales)
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis")
png("num_mQTL.png",units="in",height=5,width=6,res=1200)
ggplot(data = modified_mQTL_feature_df, aes(x = Tissue, y = num_mQTL,fill=Tissue)) +
  geom_bar(stat = "identity", width = 0.75) +
  # scale_y_continuous(trans='log10',labels = label_comma()) +
  scale_y_continuous(labels = label_comma(),expand = expansion(mult = c(0, 0.2))) +
  coord_flip() +
  labs(x = "Tissue type", y = "Number of mapped mQTLs") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = comma(num_mQTL), hjust = -0.15), size = 2.25) + 
  facet_grid(Dataset ~ ., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0),
        panel.border = element_rect(color = "black",fill = NA, size = 0.65), legend.position = "none") +
  scale_fill_manual(values = tissue_colors, guide = guide_legend(nrow = 2))
dev.off()
