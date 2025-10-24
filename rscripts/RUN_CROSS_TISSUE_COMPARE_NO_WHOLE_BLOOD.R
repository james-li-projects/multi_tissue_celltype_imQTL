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
# Step 1: Get all files with prefix "sig_" in the current working directory, excluding whole blood
files <- list.files(pattern = "^sig_.*\\.tsv$")
files <- files[!grepl("wb",files)]
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
# Step 1: Get all files with prefix "sig_" in the current working directory, excluding whole blood
files <- list.files(pattern = "^sig_.*\\.tsv$")
files <- files[!grepl("wb",files)]

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
wide_parsed_imqtl <- wide_parsed_imqtl %>% filter(tissue!="wb")


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
  GTEx_prostate_SM    = "Prostate - Smooth muscle cell"
)


##########################
# converting effect size matrix into a dataframe
beta_df<-data.frame(sig_Bhat)
# creating a correlation plot
beta_df <- na.omit(beta_df)

selected_columns <- c(
  "GTEx_colon_EC", "GTEx_colon_Epi", "GTEx_colon_Lym", "GTEx_colon_Mye", "GTEx_colon_Stromal",
  "GTEx_lung_Endo", "GTEx_lung_Epi", "GTEx_lung_Gran", "GTEx_lung_Lym", "GTEx_lung_Macro",
  "GTEx_lung_Mono", "GTEx_lung_Stromal", "GTEx_ovary_EndoC", "GTEx_ovary_IC", "GTEx_prostate_BE", "GTEx_prostate_LE", "GTEx_prostate_SM"
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
png(filename=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/all_combinations_nowb.png"),width = 8, height = 8, unit="in", res = 1200)  # High resolution
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
