library(stringr)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(Hmisc)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(Hmisc)

# set seed
set.seed(1)

#######################################
# IMPORTING Bhat ESTIMATES FOR IMQTLS #
#######################################
# importing beta hat estimates for sig imQTLs
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/b_gi")
# Step 1: Get all files with prefix "sig_" in the current working directory
files <- list.files(pattern = "^sig_.*\\.tsv$")
# Step 2: Initialize an empty list to store individual data.frames
file_list <- list()
# Step 3: Loop through each file and process it
for (file in files) {
  # Step 3a: Extract the new column name by removing specific substrings
  new_col_name <- file %>%
    str_remove("tensorQTL_imQTL_") %>%
    str_remove("^sig_") %>%
    str_remove("^sig_") %>%
    str_remove("\\.tsv$")
  # Step 3b: Read the file as a data.table
  dt <- fread(file)
  # Step 3c: Rename the columns as "pair_id" and the new column name
  colnames(dt) <- c("pair_id", new_col_name)
  # Step 3d: Add the data.table to the list
  file_list[[file]] <- dt
}
# Step 4: Join all data.tables into a single data.frame
sig_Bhat <- Reduce(function(x, y) full_join(x, y, by = "pair_id"), file_list)
sig_Bhat_df <- data.frame(sig_Bhat)
rownames(sig_Bhat) <- sig_Bhat$pair_ID
# Step 5: Convert into a matrix
sig_matrix <- as.matrix(sig_Bhat[, -"pair_id", with = FALSE])
rownames(sig_matrix) <- sig_Bhat$pair_id
sig_Bhat<-sig_matrix


#######################################
# IMPORTING Shat ESTIMATES FOR IMQTLS #
#######################################
# importing beta hat estimates for sig imQTLs
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/b_gi_se")
# Step 1: Get all files with prefix "sig_" in the current working directory
files <- list.files(pattern = "^sig_.*\\.tsv$")
# Step 2: Initialize an empty list to store individual data.frames
file_list <- list()
# Step 3: Loop through each file and process it
for (file in files) {
  # Step 3a: Extract the new column name by removing specific substrings
  new_col_name <- file %>%
    str_remove("tensorQTL_imQTL_") %>%
    str_remove("^sig_") %>%
    str_remove("^sig_") %>%
    str_remove("\\.tsv$")
  # Step 3b: Read the file as a data.table
  dt <- fread(file)
  # Step 3c: Rename the columns as "pair_id" and the new column name
  colnames(dt) <- c("pair_id", new_col_name)
  # Step 3d: Add the data.table to the list
  file_list[[file]] <- dt
}
# Step 4: Join all data.tables into a single data.frame
sig_Shat <- Reduce(function(x, y) full_join(x, y, by = "pair_id"), file_list)
sig_Shat_df <- data.frame(sig_Shat)
rownames(sig_Shat) <- sig_Shat$pair_ID
# Step 5: Convert into a matrix
sig_matrix <- as.matrix(sig_Shat[, -"pair_id", with = FALSE])
rownames(sig_matrix) <- sig_Shat$pair_id
sig_Shat<-sig_matrix


#############################################
# Identifying nominally significant effects #
#############################################
# Compute the z-scores
sig_Z <- sig_Bhat / sig_Shat

# Compute two-sided p-values
sig_Pval <- 2 * pnorm(-abs(sig_Z))

# Create a filtered version of sig_Bhat
sig_Bhat[sig_Pval >= 0.05] <- 0


###########################
# loading in reference of all imQTLs
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")



##########################
conversion_dict <- c(
  GTEx_colon_EC       = "Colon - Endothelial cell",
  GTEx_colon_Epi      = "Colon - Epithelial cell",
  GTEx_colon_Lym      = "Colon - Lymphocyte",
  GTEx_colon_Mye      = "Colon - Myeloid cell",
  GTEx_colon_Stromal  = "Colon - Stromal cell",
  GTEx_lung_Endo      = "Lung - Endothelial cell",
  GTEx_lung_Epi       = "Lung - Epithelial cell",
  GTEx_lung_Gran      = "Lung - Granulocyte",
  GTEx_lung_Lym       = "Lung - Lymphocyte",
  GTEx_lung_Macro     = "Lung - Macrophage",
  GTEx_lung_Mono      = "Lung - Monocyte",
  GTEx_lung_Stromal   = "Lung - Stromal cell",
  GTEx_ovary_EndoC    = "Ovary - Endothelial cell",
  GTEx_ovary_IC       = "Ovary - Immune cells",
  GTEx_prostate_BE    = "Prostate - Basal epithelium",
  GTEx_prostate_LE    = "Prostate - Luminal epithelium",
  GTEx_prostate_SM    = "Prostate - Smooth muscle cell",
  HEALS_wb_B          = "Whole Blood - B cell",
  HEALS_wb_CD4T       = "Whole Blood - CD4+ T cell",
  HEALS_wb_CD8T       = "Whole Blood - CD8+ T cell",
  HEALS_wb_Mono       = "Whole Blood - Monocyte",
  HEALS_wb_NK         = "Whole Blood - NK cell",
  HEALS_wb_Neutro     = "Whole Blood - Neutrophil"
)


##########################
# converting effect size matrix into a dataframe
beta_df<-data.frame(sig_Bhat)
# creating a correlation plot
beta_df <- na.omit(beta_df)

selected_columns <- c(
  "GTEx_colon_EC", "GTEx_colon_Epi", "GTEx_colon_Lym", "GTEx_colon_Mye", "GTEx_colon_Stromal",
  "GTEx_lung_Endo", "GTEx_lung_Epi", "GTEx_lung_Gran", "GTEx_lung_Lym", "GTEx_lung_Macro",
  "GTEx_lung_Mono", "GTEx_lung_Stromal", "GTEx_ovary_EndoC", "GTEx_ovary_IC", "GTEx_prostate_BE", "GTEx_prostate_LE", "GTEx_prostate_SM", "HEALS_wb_B", "HEALS_wb_CD4T", "HEALS_wb_CD8T", "HEALS_wb_Mono", "HEALS_wb_NK", "HEALS_wb_Neutro"
)

# Step 1: Filter columns
selected_columns <- intersect(colnames(beta_df), names(conversion_dict))
beta_df <- beta_df[, selected_columns]

# Step 2: Rename columns to human-readable labels
colnames(beta_df) <- conversion_dict[selected_columns]


#########################################
# Examining all imQTL variant/CpG pairs #
#########################################
# Filter relevant columns
tmp_beta_df <- beta_df #[, grep("colon|lung|ovary|prostate|blood", colnames(beta_df), ignore.case = TRUE)]

###########################
# Create a sign matrix: +1 for positive, -1 for negative, 0 for zero
sign_matrix <- sign(tmp_beta_df)
# Create an empty matrix to hold proportions
n_cols <- ncol(sign_matrix)
prop_sign_agree_mat <- matrix(NA, nrow = n_cols, ncol = n_cols)
colnames(prop_sign_agree_mat) <- colnames(sign_matrix)
rownames(prop_sign_agree_mat) <- colnames(sign_matrix)
# Compute proportion of sign agreement
for (i in 1:n_cols) {
  for (j in 1:n_cols) {
    same_sign <- sign_matrix[, i] == sign_matrix[, j]
    prop_sign_agree_mat[i, j] <- sum(same_sign, na.rm = TRUE) / nrow(sign_matrix)
  }
}
# View the matrix
round(prop_sign_agree_mat, 3)


###########################
# Load required package
library(pheatmap)
# Create sign matrix
sign_matrix <- sign(tmp_beta_df)
# Compute sign agreement proportions
n_cols <- ncol(sign_matrix)
prop_sign_agree_mat <- matrix(NA, nrow = n_cols, ncol = n_cols)
colnames(prop_sign_agree_mat) <- colnames(sign_matrix)
rownames(prop_sign_agree_mat) <- colnames(sign_matrix)
for (i in 1:n_cols) {
  for (j in 1:n_cols) {
    same_sign <- sign_matrix[, i] == sign_matrix[, j]
    prop_sign_agree_mat[i, j] <- sum(same_sign, na.rm = TRUE) / nrow(sign_matrix)
  }
}
# Optional: Round for readability
rounded_mat <- round(prop_sign_agree_mat, 2)

# Save high-quality PNG heatmap
png(filename=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/all_combinations.png"),width = 8, height = 8, unit="in", res = 1200)  # High resolution
# Plot
pheatmap(rounded_mat,
         color = colorRampPalette(c("white", "#FFFF99","#FFFF80","orangered","red"))(100),
         border_color = NA,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 6,
         fontsize_row = 8,
         fontsize_col = 8,
         clustering_method = "complete",
         main = "Proportion of Shared Effects")
dev.off()


#####################################
# Examining subpopulations of cells #
#####################################
# cell type lists
list_lymphoid <- c("B", "CD4T", "CD8T", "Lym", "NK")
list_myeloid <- c("Gran", "Macro", "Mono", "Mye", "Neutro")
list_all_immune_cell <- c("B", "CD4T", "CD8T", "Lym", "NK", "Gran", "Macro", "Mono", "Mye", "Neutro", "IC")
list_epithelial <- c("BE", "Epi", "LE")
# combine lists
cell_type_lists <- list(
  lymphoid = list_lymphoid,
  myeloid = list_myeloid,
  all_immune_cell = list_all_immune_cell,
  epithelial = list_epithelial
)

for (list_name in names(cell_type_lists)) {
  filtered_imqtl <- wide_parsed_imqtl %>%
    filter(celltype %in% cell_type_lists[[list_name]]) %>%
    mutate(pair_id = paste(phenotype_id, variant_id, sep = "_"))
  
  ###########################
  # Filter relevant columns
  tmp_beta_df <- beta_df[rownames(beta_df) %in% filtered_imqtl$pair_id,]
  tmp_beta_df <- tmp_beta_df
  
  ###########################
  # Create a sign matrix: +1 for positive, -1 for negative, 0 for zero
  sign_matrix <- sign(tmp_beta_df)
  # Create an empty matrix to hold proportions
  n_cols <- ncol(sign_matrix)
  prop_sign_agree_mat <- matrix(NA, nrow = n_cols, ncol = n_cols)
  colnames(prop_sign_agree_mat) <- colnames(sign_matrix)
  rownames(prop_sign_agree_mat) <- colnames(sign_matrix)
  # Compute proportion of sign agreement
  for (i in 1:n_cols) {
    for (j in 1:n_cols) {
      same_sign <- sign_matrix[, i] == sign_matrix[, j]
      prop_sign_agree_mat[i, j] <- sum(same_sign, na.rm = TRUE) / nrow(sign_matrix)
    }
  }
  # View the matrix
  round(prop_sign_agree_mat, 3)
  
  
  ###########################
  # Load required package
  library(pheatmap)
  # Create sign matrix
  sign_matrix <- sign(tmp_beta_df)
  # Compute sign agreement proportions
  n_cols <- ncol(sign_matrix)
  prop_sign_agree_mat <- matrix(NA, nrow = n_cols, ncol = n_cols)
  colnames(prop_sign_agree_mat) <- colnames(sign_matrix)
  rownames(prop_sign_agree_mat) <- colnames(sign_matrix)
  for (i in 1:n_cols) {
    for (j in 1:n_cols) {
      same_sign <- sign_matrix[, i] == sign_matrix[, j]
      prop_sign_agree_mat[i, j] <- sum(same_sign, na.rm = TRUE) / nrow(sign_matrix)
    }
  }
  # Optional: Round for readability
  rounded_mat <- round(prop_sign_agree_mat, 2)
  
  # Save high-quality PNG heatmap
  png(filename=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/all_combinations_",list_name,".png"),width = 8, height = 8, unit="in", res = 1200)  # High resolution
  # Plot
  pheatmap(rounded_mat,
           color = colorRampPalette(c("white", "#FFFF99","#FFFF80","orangered","red"))(100),
           border_color = NA,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 6,
           fontsize_row = 8,
           fontsize_col = 8,
           clustering_method = "complete",
           main = "Proportion of Shared Effects")
  dev.off()    
}
































































































# combination list
combination_list_above20 <- as.vector((data.frame(table(wide_parsed_imqtl$combination)) %>% filter(Freq>100))$Var1)

for (current_filter_combination in combination_list_above20) {
  filtered_imqtl <- wide_parsed_imqtl %>%
    filter(combination==current_filter_combination) %>%
    mutate(pair_id = paste(phenotype_id, variant_id, sep = "_"))
  
  
  # Load necessary libraries
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(Hmisc)  # For rcorr function
  
  # Filter relevant columns
  tmp_beta_df <- beta_df[rownames(beta_df) %in% filtered_imqtl$pair_id,]
  tmp_beta_df <- tmp_beta_df[, grep("colon|lung|ovary|prostate|HEALS_wb", colnames(beta_df), ignore.case = TRUE)]

  
  # Compute Spearman correlation matrix and p-values
  rcorr_result <- Hmisc::rcorr(as.matrix(tmp_beta_df), type = "spearman")
  cor_matrix <- rcorr_result$r
  pval_matrix <- rcorr_result$P
  
  # Create matrix of asterisks: * for p < 0.05, ** for p < 0.05/253
  star_matrix <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
  rownames(star_matrix) <- rownames(pval_matrix)
  colnames(star_matrix) <- colnames(pval_matrix)
  
  # Fill in asterisks based on significance thresholds
  for (i in 1:nrow(pval_matrix)) {
    for (j in 1:ncol(pval_matrix)) {
      if (i != j) {
        if (pval_matrix[i, j] < 0.05 / 253) {
          star_matrix[i, j] <- "*"
        } else if (pval_matrix[i, j] < 0.05) {
          star_matrix[i, j] <- ""
        }
      }
    }
  }
  
  # Save as high quality PNG
  png(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/all_combinations_",current_filter_combination,".png"), 
      width = 2400, height = 2000, res = 300)
  
  # Plot the heatmap with asterisks in cells
  pheatmap(cor_matrix,
           display_numbers = star_matrix,
           main = "Correlation Matrix of Beta Values Across Tissues",
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 8,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           border_color = NA)
  
  # Close the graphics device
  dev.off()
  
}















































###############################
# ASSEMBLING MASHR INPUT LIST #
###############################
# importing significant vs random imQTL pairs into mashR
strong.subset = mash_set_data(sig_Bhat,sig_Shat) 
random.subset = mash_set_data(random_Bhat,random_Shat) 

# calculating null correlation matrix
data.temp=random.subset
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

# finalizing data that is tuned based on this null correlation matrix
data.random = mash_set_data(random_Bhat,random_Shat,V=Vhat)
data.strong = mash_set_data(sig_Bhat,sig_Shat, V=Vhat)

# set up data-driven covariances
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)


#################
# RUNNING MASHR #
#################
# fit mash to the random subset using the data-driven covariances
Sys.time()
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
Sys.time()

# performing cross_tissue_compare to compute posterior summaries
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
Sys.time()
cross_tissue_compare_pairwise_tissue_sharing<-cross_tissue_compare::get_pairwise_sharing(m2)
Sys.time()


#####################
# SAVE OUTPUT FILES #
#####################
save(data.random,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","data.random",".RData"))
save(data.strong,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","data.strong",".RData"))
save(U.pca,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","U.pca",".RData"))
save(U.ed,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","U.ed",".RData"))
save(U.c,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","U.c",".RData"))
save(m,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","m",".RData"))
save(m2,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","m2",".RData"))
save(cross_tissue_compare_pairwise_tissue_sharing,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","cross_tissue_compare_pairwise_tissue_sharing",".RData"))

#save(,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/","",".RData"))
