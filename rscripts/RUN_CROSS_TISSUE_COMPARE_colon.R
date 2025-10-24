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

# filtering for tissue of interest in this DF and sig_Bhat
wide_parsed_imqtl <- wide_parsed_imqtl %>% filter(tissue=="colon") %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_"),extract_col_id=paste0("GTEx_",tissue,"_",celltype))
sig_Bhat <- sig_Bhat[wide_parsed_imqtl$pair_id,unique(wide_parsed_imqtl$extract_col_id)]

##########################
conversion_dict <- c(
  GTEx_colon_EC       = "Colon - Endothelial cell",
  GTEx_colon_Epi      = "Colon - Epithelial cell",
  GTEx_colon_Lym      = "Colon - Lymphocyte",
  GTEx_colon_Mye      = "Colon - Myeloid cell",
  GTEx_colon_Stromal  = "Colon - Stromal cell"
)


##########################
# converting effect size matrix into a dataframe
beta_df<-data.frame(sig_Bhat)
# creating a correlation plot
beta_df <- na.omit(beta_df)

selected_columns <- c(
  "GTEx_colon_EC", "GTEx_colon_Epi", "GTEx_colon_Lym", "GTEx_colon_Mye", "GTEx_colon_Stromal"
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
tmp_beta_df <- beta_df

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

# Only retain cell types in names
colnames(rounded_mat) <- gsub("Colon - ","",colnames(rounded_mat))
rownames(rounded_mat) <- gsub("Colon - ","",rownames(rounded_mat))

# Save high-quality PNG heatmap
png(filename=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/prop_shared_effects_colon.png"),width = 5, height = 5, unit="in", res = 1200)  # High resolution
# Plot
pheatmap(rounded_mat,
         color = colorRampPalette(c("white", "#FFFF99","#FFFF80","orangered","red"))(100),
         border_color = NA,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 13,
         fontsize_row = 13,
         fontsize_col = 13,
         clustering_method = "complete",
         angle_col = 45)
dev.off()
