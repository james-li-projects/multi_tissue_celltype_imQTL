library(ashr)
library(cross_tissue_compare)
library(stringr)
library(data.table)
library(dplyr)
library(stringr)

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
# importing beta hat estimates for random imQTLs
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/b_gi")
# Step 1: Get all files with prefix "random_" in the current working directory
files <- list.files(pattern = "^random_.*\\.tsv$")
# Step 2: Initialize an empty list to store individual data.frames
file_list <- list()
# Step 3: Loop through each file and process it
for (file in files) {
  # Step 3a: Extract the new column name by removing specific substrings
  new_col_name <- file %>%
    str_remove("tensorQTL_imQTL_") %>%
    str_remove("^random_") %>%
    str_remove("^random_") %>%
    str_remove("\\.tsv$")
  # Step 3b: Read the file as a data.table
  dt <- fread(file)
  # Step 3c: Rename the columns as "pair_id" and the new column name
  colnames(dt) <- c("pair_id", new_col_name)
  # Step 3d: Add the data.table to the list
  file_list[[file]] <- dt
}
# Step 4: Join all data.tables into a single data.frame
random_Bhat <- Reduce(function(x, y) full_join(x, y, by = "pair_id"), file_list)
random_Bhat_df <- data.frame(random_Bhat)
rownames(random_Bhat) <- random_Bhat$pair_ID
# Step 5: Convert into a matrix
random_matrix <- as.matrix(random_Bhat[, -"pair_id", with = FALSE])
rownames(random_matrix) <- random_Bhat$pair_id
random_Bhat<-random_matrix


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

#######################################
# importing beta hat estimates for random imQTLs
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/b_gi_se")
# Step 1: Get all files with prefix "random_" in the current working directory
files <- list.files(pattern = "^random_.*\\.tsv$")
# Step 2: Initialize an empty list to store individual data.frames
file_list <- list()
# Step 3: Loop through each file and process it
for (file in files) {
  # Step 3a: Extract the new column name by removing specific substrings
  new_col_name <- file %>%
    str_remove("tensorQTL_imQTL_") %>%
    str_remove("^random_") %>%
    str_remove("^random_") %>%
    str_remove("\\.tsv$")
  # Step 3b: Read the file as a data.table
  dt <- fread(file)
  # Step 3c: Rename the columns as "pair_id" and the new column name
  colnames(dt) <- c("pair_id", new_col_name)
  # Step 3d: Add the data.table to the list
  file_list[[file]] <- dt
}
# Step 4: Join all data.tables into a single data.frame
random_Shat <- Reduce(function(x, y) full_join(x, y, by = "pair_id"), file_list)
random_Shat_df <- data.frame(random_Shat)
rownames(random_Shat) <- random_Shat$pair_ID
# Step 5: Convert into a matrix
random_matrix <- as.matrix(random_Shat[, -"pair_id", with = FALSE])
rownames(random_matrix) <- random_Shat$pair_id
random_Shat<-random_matrix


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
