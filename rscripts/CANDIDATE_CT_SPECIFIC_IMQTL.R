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

# identifying candidate cell-type specific mQTLs
CT_specific_imQTL_df <- data.frame()
for (tmp_tissue in c("wb","colon","lung")) {
  # setting dataset based on tissue analyzed, since wb imQTLs are uniquely identified in HEALS
  tmp_dataset <- ifelse(tmp_tissue=="wb","HEALS","GTEx")
  
  # identifying all columns for cell types in the current tissue analyzed 
  tissue_wide_parsed_imqtl <- wide_parsed_imqtl %>% filter(tissue==tmp_tissue) %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_"),extract_col_id=paste0(tmp_dataset,"_",tissue,"_",celltype))
  tissue_specific_celltype_col_id <- unique(tissue_wide_parsed_imqtl$extract_col_id)
  
  # iterating through each cell type to identify candidate cell type specific imQTLs
  for (focus_combination in tissue_specific_celltype_col_id) {
    print(focus_combination)
    
    # focusing our effect size matrix on those variant-cpg pairs only identified as imQTLs for the current cell type 
    tmp_combination_pair_id <- (tissue_wide_parsed_imqtl %>% filter(extract_col_id==focus_combination))$pair_id
    tmp_sig_Bhat <- sig_Bhat[tmp_combination_pair_id,tissue_specific_celltype_col_id, drop = FALSE]
    
    # converting effect size matrix into a dataframe
    beta_df <- as.data.frame(tmp_sig_Bhat)
    # removing NA values
    beta_df <- na.omit(beta_df)
    # filter columns for those cell types in the current tissue
    beta_df <- beta_df[, tissue_specific_celltype_col_id]
    
    # identify candidate cell-type specific mQTLs
    # Get sign of each value
    sign_mat <- sign(beta_df)
    
    # Separate target column and others
    target_vals <- sign_mat[, focus_combination]
    other_vals <- sign_mat[, setdiff(colnames(sign_mat), focus_combination)]
    
    # Apply generalized filtering
    keep <- (target_vals ==  1 & apply(other_vals, 1, function(x) all(x %in% c(-1, 0)))) |
      (target_vals == -1 & apply(other_vals, 1, function(x) all(x %in% c(1, 0))))
    
    filtered_sign_mat <- sign_mat[keep, ]
    print(nrow(filtered_sign_mat))
    print(nrow(beta_df))
    
    # store imQTL results for cell-type specific hits
    tmp_CT_specific_imQTL_df <- tissue_wide_parsed_imqtl %>% filter(pair_id %in% rownames(filtered_sign_mat),extract_col_id==focus_combination)
    CT_specific_imQTL_df <- rbind(CT_specific_imQTL_df,tmp_CT_specific_imQTL_df)
  }
}



