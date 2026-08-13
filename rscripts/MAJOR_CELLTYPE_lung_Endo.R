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
library(ggnewscale)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

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

###########################
# loading in reference of all imQTLs
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")

# identifying columns to cell types to extract
extract_col_id=paste0("GTEx_",unique((wide_parsed_imqtl %>% filter(tissue=="lung"))$combination))

# filtering for tissue of interest in this DF and sig_Bhat
wide_parsed_imqtl <- wide_parsed_imqtl %>% filter(combination=="lung_Endo") %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_"),extract_col_id=paste0("GTEx_",tissue,"_",celltype))
sig_Bhat <- sig_Bhat[wide_parsed_imqtl$pair_id,extract_col_id]

##########################
conversion_dict <- c(
  GTEx_lung_Endo      = "Lung - Endothelial cell",
  GTEx_lung_Epi       = "Lung - Epithelial cell",
  GTEx_lung_Gran      = "Lung - Granulocyte",
  GTEx_lung_Lym       = "Lung - Lymphocyte",
  GTEx_lung_Macro     = "Lung - Macrophage",
  GTEx_lung_Mono      = "Lung - Monocyte",
  GTEx_lung_Stromal   = "Lung - Stromal cell"
)

##########################
# converting effect size matrix into a dataframe
beta_df<-data.frame(sig_Bhat)
# creating a correlation plot
beta_df <- na.omit(beta_df)

selected_columns <- c(
  "GTEx_lung_Endo", "GTEx_lung_Epi", "GTEx_lung_Gran", "GTEx_lung_Lym", "GTEx_lung_Macro",
  "GTEx_lung_Mono", "GTEx_lung_Stromal"
)

# Step 1: Filter columns
selected_columns <- intersect(colnames(beta_df), names(conversion_dict))
beta_df <- beta_df[, selected_columns]

# Step 2: Rename columns to human-readable labels
colnames(beta_df) <- conversion_dict[selected_columns]

#########################################################
# Examine effect sizes in comparison to major cell type
major_celltype_beta_df <- beta_df
flip_rows <- major_celltype_beta_df[["Lung - Endothelial cell"]] < 0
major_celltype_beta_df[flip_rows, ] <- -beta_df[flip_rows, ]
major_celltype_beta_df[["Lung - Endothelial cell"]] <- abs(major_celltype_beta_df[["Lung - Endothelial cell"]])

# Reshape to long format and clean cell type names
long_df <- major_celltype_beta_df %>%
  tibble::rownames_to_column("cpg_id") %>%
  pivot_longer(cols = -cpg_id, names_to = "cell_type", values_to = "effect_size") %>%
  mutate(cell_type = gsub("Lung - ", "", cell_type))

###################################
# creating violin plots
library(ggplot2)
library(dplyr)
library(ggnewscale)

# Specify your target cell type here
target_celltype <- "Endothelial cell"

# Ensure cell_type is a factor and get levels with target first
all_celltypes <- unique(long_df$cell_type)
cell_type_order <- c(target_celltype, setdiff(all_celltypes, target_celltype))
long_df <- long_df %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_order))

# Define violin colors: gold for target, navy for others
violin_colors <- setNames(
  ifelse(cell_type_order == target_celltype, "gold", "navy"),
  cell_type_order
)

# Boxplot colors: all white
box_colors <- setNames(rep("white", length(cell_type_order)), cell_type_order)

# Plot
p <- ggplot(long_df, aes(x = cell_type, y = effect_size)) +
  geom_violin(aes(fill = cell_type), trim = FALSE, color = "black", alpha = 0.7) +
  scale_fill_manual(values = violin_colors) +
  guides(fill = "none") +
  new_scale("fill") +
  geom_boxplot(aes(fill = cell_type), width = 0.12, outlier.shape = NA, alpha = 1, color = "black") +
  scale_fill_manual(values = box_colors) +
  theme_classic(base_size = 16) +
  labs(
    x = "",
    y = "Interaction Term Estimate"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6)

# Save plot
ggsave(
  filename = paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cross_tissue_compare_output/violinplot_", target_celltype, ".png"),
  plot = p,
  width = 4.5,
  height = 5,
  dpi = 600
)
