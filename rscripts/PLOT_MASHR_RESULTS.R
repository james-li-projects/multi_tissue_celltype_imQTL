library(tidyr)

imQTL_feature_df <- data.frame()

for (i in 1:ncol(m2$result$lfsr)) {
  lfsr_0.05_indices <- (which((m2$result$lfsr[,i])<0.05))
  sig_pairs_lfsr_0.05<-rownames(m2$result$lfsr)[lfsr_0.05_indices]
  sig_cpgs_lfsr_0.05 <- sapply(strsplit(sig_pairs_lfsr_0.05, "_"), `[`, 1)
  current_combination_name=colnames(m2$result$lfsr)[i]
  print(paste(
    current_combination_name,
    length(unique(sig_cpgs_lfsr_0.05)))
  )
  
  imQTL_feature_df <- rbind(
    imQTL_feature_df,
    data.frame(combination=current_combination_name,num_imQTL=length(unique(sig_cpgs_lfsr_0.05)))                        
  )
}
imQTL_feature_df <- imQTL_feature_df %>% separate(combination,sep="_",into=c("Dataset","Tissue","Celltype"),remove=F)


#######################################
# making the tissue names in the imQTL_feature_df more beautiful
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="breast"] <- "Breast"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="colon"] <- "Colon"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="lung"] <- "Lung"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="kidney"] <- "Kidney"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="prostate"] <- "Prostate"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="wb"] <- "Whole Blood"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="ovary"] <- "Ovary"
# making the cell type names in the imQTL_feature_df more beautiful
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Basal"] <- "Basal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Luminal"] <- "Luminal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Epi"] <- "Epithelial cell"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="BE"] <- "Basal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="LE"] <- "Luminal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="SM"] <- "Smooth muscle cell"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="EC"] <- "Endothelial cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Fat"] <- "Adipocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Fib"] <- "Fibroblast"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Lym"] <- "Lymphocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="MP"] <- "Macrophage"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Macro"] <- "Macrophage"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Mono"] <- "Monocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Mye"] <- "Myeloid cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Gran"] <- "Granulocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Stromal"] <- "Stromal cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Endo"] <- "Endothelial cell"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="B"] <- "B cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="CD4T"] <- "CD4+ T Cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="CD8T"] <- "CD8+ T Cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="NK"] <- "NK cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Neutro"] <- "Neutrophil"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="EndoC"] <- "Endothelial cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="IC"] <- "Immune cells"

imQTL_feature_df$imQTL_vec <- paste(imQTL_feature_df$Tissue,imQTL_feature_df$Celltype,sep=" - ")


#######################################
# generating plots
library(ggplot2)
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/mashr/mashr_output")
png("num_mashr_imQTL.png",units="in",height=5,width=6,res=1200)
ggplot(data = imQTL_feature_df, aes(x = imQTL_vec, y = num_imQTL,fill=Tissue)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_y_continuous(trans='log10') +
  coord_flip() +
  labs(x = "Cell Type", y = "Number of mapped imQTLs") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = num_imQTL, hjust = 1.1), size = 2.25) + 
  #facet_wrap(~Dataset,ncol=1,shrink=T) 
  facet_grid(Dataset ~ ., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0),
        panel.border = element_rect(color = "black",fill = NA, size = 0.65))
dev.off()






#############################################
# List of groups of names to create smaller correlation matrices
selected_groups <- list(
  group1 = colnames(m2$result$lfsr)[grepl("breast",colnames(m2$result$lfsr))],
  group2 = colnames(m2$result$lfsr)[grepl("kidney",colnames(m2$result$lfsr))],
  group3 = colnames(m2$result$lfsr)[grepl("colon",colnames(m2$result$lfsr))],
  group4 = colnames(m2$result$lfsr)[grepl("lung",colnames(m2$result$lfsr))],
  group5 = colnames(m2$result$lfsr)[grepl("ovary",colnames(m2$result$lfsr))],
  group6 = colnames(m2$result$lfsr)[grepl("prostate",colnames(m2$result$lfsr))],
  group7 = colnames(m2$result$lfsr)[grepl("HEALS_wb",colnames(m2$result$lfsr))],
  group8 = colnames(m2$result$lfsr)[grepl("GTEx_wb",colnames(m2$result$lfsr))]
)
# Function to create smaller correlation matrices
create_sub_matrices <- function(corr_matrix, groups) {
  lapply(groups, function(group) corr_matrix[group, group, drop = FALSE])
}

# Create smaller correlation matrices
small_corr_matrices <- create_sub_matrices(mashr_pairwise_tissue_sharing, selected_groups)


######################
# Load required library
library(corrplot)

# Define the output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/mashr/mashr_output"

######################
# Define the file name for the corr plot without mashr 
file_name <- file.path(output_dir, "mashr_no_mashr_plot.png")

# Calculate correlation matrix while handling NAs
cor_matrix <- cor(sig_Bhat, use = "pairwise.complete.obs")

# Save the plot as a high-resolution PNG
png(file_name, width = 5000, height = 5000, res = 300) # Adjust resolution as needed
corrplot(
  cor_matrix, 
  method = "color", 
  type = "upper", 
  tl.col = "black", 
  tl.srt = 45
)
dev.off()









# plotting submatrices for GTEx tissues
tissue_list=c("colon","lung","ovary","prostate")
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/imQTL/top_assoc")
imqtl_file_list = list.files()
imqtl_file_list = imqtl_file_list[grepl(".cis_qtl_top_assoc.txt.gz",imqtl_file_list)]
for (current_tissue in tissue_list) {
  print(current_tissue)
  tmp_pair_list <- c()
  tmp_imqtl_file_list <- imqtl_file_list[grepl(current_tissue,imqtl_file_list)]
  for (tmp_imqtl_file in tmp_imqtl_file_list) {
    tmp_imqtl_df <- fread(tmp_imqtl_file) %>% filter(pval_adj_bh<0.05) %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_"))
    tmp_pair_list <- c(tmp_pair_list,tmp_imqtl_df$pair_id)
  }
  
  # filtering matrix for pair IDs and tissue type as well as not for tissue type
  tissue_filtered_sig_Bhat <- sig_Bhat[rownames(sig_Bhat) %in% tmp_pair_list, grepl(current_tissue,colnames(sig_Bhat))]
  all_filtered_sig_Bhat <- sig_Bhat[rownames(sig_Bhat) %in% tmp_pair_list, ]
  
  # Calculate correlation matrix while handling NAs
  tissue_cor_matrix <- cor(tissue_filtered_sig_Bhat, use = "pairwise.complete.obs")
  all_cor_matrix <- cor(all_filtered_sig_Bhat, use = "pairwise.complete.obs")

  # plot corrplot for a specific tissue
  png(paste0(output_dir,"/tissue_corrplot_",current_tissue,".png"), width = 5000, height = 5000, res = 300) # Adjust resolution as needed
  corrplot(
    tissue_cor_matrix, 
    method = "color", 
    type = "upper", 
    tl.col = "black", 
    tl.srt = 45,
    tl.cex = 2,      # Increase text label size
    number.cex = 1.5 # Increase number size 
  )
  dev.off()
  
  # plot corrplot across all tissues
  png(paste0(output_dir,"/all_corrplot_",current_tissue,".png"), width = 5000, height = 5000, res = 300) # Adjust resolution as needed
  corrplot(
    all_cor_matrix, 
    method = "color", 
    type = "upper", 
    tl.col = "black", 
    tl.srt = 45
  )
  dev.off()
}



# plotting submatrices for HEALS whole blood
tissue_list=c("wb")
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc")
imqtl_file_list = list.files()
imqtl_file_list = imqtl_file_list[grepl(".cis_qtl_top_assoc.txt.gz",imqtl_file_list)]
for (current_tissue in tissue_list) {
  print(current_tissue)
  tmp_pair_list <- c()
  tmp_imqtl_file_list <- imqtl_file_list[grepl(current_tissue,imqtl_file_list)]
  for (tmp_imqtl_file in tmp_imqtl_file_list) {
    tmp_imqtl_df <- fread(tmp_imqtl_file) %>% filter(pval_adj_bh<0.05) %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_"))
    tmp_pair_list <- c(tmp_pair_list,tmp_imqtl_df$pair_id)
  }
  
  # filtering matrix for pair IDs and tissue type as well as not for tissue type
  tissue_filtered_sig_Bhat <- sig_Bhat[rownames(sig_Bhat) %in% tmp_pair_list, grepl(current_tissue,colnames(sig_Bhat))]
  all_filtered_sig_Bhat <- sig_Bhat[rownames(sig_Bhat) %in% tmp_pair_list, ]
  
  # Calculate correlation matrix while handling NAs
  tissue_cor_matrix <- cor(tissue_filtered_sig_Bhat, use = "pairwise.complete.obs")
  all_cor_matrix <- cor(all_filtered_sig_Bhat, use = "pairwise.complete.obs")
  
  # plot corrplot for a specific tissue
  png(paste0(output_dir,"/tissue_corrplot_",current_tissue,".png"), width = 5000, height = 5000, res = 300) # Adjust resolution as needed
  corrplot(
    tissue_cor_matrix, 
    method = "color", 
    type = "upper", 
    tl.col = "black", 
    tl.srt = 45,
    tl.cex = 2,      # Increase text label size
    number.cex = 1.5 # Increase number size 
  )
  dev.off()
  
  # plot corrplot across all tissues
  png(paste0(output_dir,"/all_corrplot_",current_tissue,".png"), width = 5000, height = 5000, res = 300) # Adjust resolution as needed
  corrplot(
    all_cor_matrix, 
    method = "color", 
    type = "upper", 
    tl.col = "black", 
    tl.srt = 45
  )
  dev.off()
}























######################
# Define the file name for the Vhat plot
file_name <- file.path(output_dir, "mashr_Vhat_plot.png")

# Save the plot as a high-resolution PNG
png(file_name, width = 5000, height = 5000, res = 300) # Adjust resolution as needed
corrplot(
  Vhat, 
  method = "color", 
  type = "upper", 
  tl.col = "black", 
  tl.srt = 45
)
dev.off()


######################
# Define the file name for the mashr_pairwise_tissue_sharing_plot plot
file_name <- file.path(output_dir, "mashr_pairwise_tissue_sharing_plot.png")

# Save the plot as a high-resolution PNG
png(file_name, width = 5000, height = 5000, res = 300) # Adjust resolution as needed
corrplot(
  mashr_pairwise_tissue_sharing, 
  method = "color", 
  type = "upper", 
  tl.col = "black", 
  tl.srt = 45
)
dev.off()

######################
# Iterate over each matrix in the list and create a plot
for (i in seq_along(small_corr_matrices)) {
  # Extract the current matrix
  matrix <- small_corr_matrices[[i]]
  
  # Define the file name
  file_name <- file.path(output_dir, paste0("correlation_plot_", i, ".png"))
  
  # Save the plot as a high-resolution PNG
  png(file_name, width = 3000, height = 3000, res = 300) # Adjust resolution as needed
  corrplot(
    matrix, 
    method = "color", 
    type = "upper", 
    tl.col = "black", 
    tl.srt = 45, 
    addCoef.col = "black", # Add correlation coefficients as numbers
    number.cex = 1.2       # Increase the size of the numbers
  )
  dev.off()
}


