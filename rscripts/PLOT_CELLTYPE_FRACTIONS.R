library(tidyverse)

long_estF_list <- list()

# ---- Step 1: GTEx (excluding wb) ----
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac")
gtex_files <- list.files(pattern = "^estF\\..*\\.RData$")
gtex_files <- gtex_files[!grepl("estF\\.wb\\.RData$", gtex_files)]

for (file in gtex_files) {
  load(file)
  tissue <- str_extract(file, "(?<=estF\\.)[^.]+(?=\\.RData)")
  var <- ls(pattern = "^estF\\.")
  df <- get(var)
  long_estF_list[[paste0("GTEx_", tissue)]] <- df %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(-sample_id, names_to = "celltype", values_to = "Value") %>%
    mutate(tissue = tissue, dataset = "GTEx")
  rm(list = var)
}

# ---- Step 2: HEALS (include all) ----
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac")
heals_files <- list.files(pattern = "^estF\\..*\\.RData$")

for (file in heals_files) {
  load(file)
  tissue <- str_extract(file, "(?<=estF\\.)[^.]+(?=\\.RData)")
  var <- ls(pattern = "^estF\\.")
  df <- get(var)
  long_estF_list[[paste0("HEALS_", tissue)]] <- df %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(-sample_id, names_to = "celltype", values_to = "Value") %>%
    mutate(tissue = tissue, dataset = "HEALS")
  rm(list = var)
}

# assembling DF for all cell type fractions from all samples from all tissues
all_ct_df <- bind_rows(long_estF_list)

# making celltype and tissue names more beautiful
for (j in 1:nrow(all_ct_df)) {
  tmp_ct_df <- all_ct_df[j,]
  
  # Extract unprocessed names
  unprocessed_tissue <- as.character(all_ct_df[j,"tissue"])
  unprocessed_celltype <- as.character(all_ct_df[j,"celltype"])
  
  # Create lookup vectors
  tissue_map <- c(
    "breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
    "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary"
  )
  celltype_map <- c(
    "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
    "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
    "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
    "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
    "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
    "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
    "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
  )
  # Map to beautiful names
  processed_tissue <- tissue_map[[unprocessed_tissue]]
  processed_celltype <- celltype_map[[unprocessed_celltype]]
  all_ct_df[j,"tissue"] <- processed_tissue
  all_ct_df[j,"celltype"] <- processed_celltype
}

# add HEALS label to whole blood
all_ct_df <- all_ct_df %>%
  mutate(tissue=ifelse(tissue=="Whole Blood","Whole Blood (HEALS)",tissue))

##########################
library(ggplot2)
library(dplyr)
library(forcats)

# Custom tissue colors
tissue_colors <- c(
  "Breast" = "turquoise3",
  "Lung" = "yellowgreen",
  "Colon" = "sienna",
  "Ovary" = "pink3",
  "Prostate" = "lightgray",
  "Whole Blood (HEALS)" = "magenta3",
  "Kidney" = "turquoise2"
)

# Prepare data
all_ct_df <- all_ct_df %>%
  mutate(
    tissue = factor(tissue, levels = names(tissue_colors)),
    combo = paste(celltype, tissue, sep = " - "),
    combo = factor(combo, levels = unique(combo)),
    celltype_display = celltype
  )

# File path
output_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/all_celltype_boxplots.png"

# Save PNG
png(filename = output_path, width = 3800, height = 1000, res = 300)

# Plot
ggplot(all_ct_df, aes(x = combo, y = Value, fill = tissue)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, linewidth = 0.25, color = "black") +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  scale_x_discrete(labels = setNames(all_ct_df$celltype_display, all_ct_df$combo)) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "Cell Type Proportion",
    fill = "Tissue"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    strip.background = element_rect(color = "black", fill = "grey90", linewidth = 0.6),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.spacing.x = unit(0.6, "lines"),
    legend.position = "none",
    legend.box = "horizontal"
  ) #+
  #guides(fill = guide_legend(title.position = "top", nrow = 1))

dev.off()


###############################################
###############################################
###############################################
###############################################
###############################################
# plotting correlation matrices of cell types
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(tidyverse)

long_estF_list <- list()

# ---- Step 1: GTEx (excluding wb) ----
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac")
gtex_files <- list.files(pattern = "^estF\\..*\\.RData$")
gtex_files <- gtex_files[!grepl("estF\\.wb\\.RData$", gtex_files)]

for (file in gtex_files) {
  load(file)
  tissue <- str_extract(file, "(?<=estF\\.)[^.]+(?=\\.RData)")
  var <- ls(pattern = "^estF\\.")
  df <- get(var)
  long_estF_list[[paste0("GTEx_", tissue)]] <- df %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(-sample_id, names_to = "celltype", values_to = "Value") %>%
    mutate(tissue = tissue, dataset = "GTEx")
  rm(list = var)
}

# ---- NEW: Add GTEx Whole Blood ----
if (file.exists("estF.wb.RData")) {
  load("estF.wb.RData")
  df <- get(ls(pattern = "^estF\\."))
  long_estF_list[["GTEx_wb"]] <- df %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(-sample_id, names_to = "celltype", values_to = "Value") %>%
    mutate(tissue = "wb", dataset = "GTEx")
  rm(list = ls(pattern = "^estF\\."))
}

# ---- Step 2: HEALS (include all) ----
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac")
heals_files <- list.files(pattern = "^estF\\..*\\.RData$")

for (file in heals_files) {
  load(file)
  tissue <- str_extract(file, "(?<=estF\\.)[^.]+(?=\\.RData)")
  var <- ls(pattern = "^estF\\.")
  df <- get(var)
  long_estF_list[[paste0("HEALS_", tissue)]] <- df %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(-sample_id, names_to = "celltype", values_to = "Value") %>%
    mutate(tissue = tissue, dataset = "HEALS")
  rm(list = var)
}

# ---- Combine all data ----
all_ct_df <- bind_rows(long_estF_list)

# ---- Make celltype and tissue names more readable ----
for (j in 1:nrow(all_ct_df)) {
  tmp_ct_df <- all_ct_df[j,]
  
  unprocessed_tissue <- as.character(all_ct_df[j,"tissue"])
  unprocessed_celltype <- as.character(all_ct_df[j,"celltype"])
  
  tissue_map <- c(
    "breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
    "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary"
  )
  celltype_map <- c(
    "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
    "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
    "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
    "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
    "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
    "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
    "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
  )
  
  processed_tissue <- tissue_map[[unprocessed_tissue]]
  processed_celltype <- celltype_map[[unprocessed_celltype]]
  all_ct_df[j,"tissue"] <- processed_tissue
  all_ct_df[j,"celltype"] <- processed_celltype
}

# ---- Relabel Whole Blood by dataset ----
all_ct_df <- all_ct_df %>%
  mutate(tissue = case_when(
    tissue == "Whole Blood" & dataset == "HEALS" ~ "Whole Blood (HEALS)",
    tissue == "Whole Blood" & dataset == "GTEx"  ~ "Whole Blood (GTEx)",
    TRUE ~ tissue
  ))

# ---- Custom tissue colors ----
tissue_colors <- c(
  "Breast" = "turquoise3",
  "Lung" = "yellowgreen",
  "Colon" = "sienna",
  "Ovary" = "pink3",
  "Prostate" = "lightgray",
  "Whole Blood (HEALS)" = "magenta3",
  "Whole Blood (GTEx)" = "magenta3",
  "Kidney" = "turquoise2"
)

# ---- Final formatting ----
all_ct_df <- all_ct_df %>%
  mutate(
    tissue = factor(tissue, levels = names(tissue_colors)),
    combo = paste(celltype, tissue, sep = " - "),
    combo = factor(combo, levels = unique(combo)),
    celltype_display = celltype
  )

# Output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/"

# Loop through each tissue
for (t in unique(all_ct_df$tissue)) {
  
  # Prepare data
  tissue_df <- all_ct_df %>%
    filter(tissue == t) %>%
    select(sample_id, celltype_display, Value) %>%
    pivot_wider(names_from = celltype_display, values_from = Value) %>%
    column_to_rownames("sample_id")
  
  # Drop columns with all NAs
  tissue_df <- tissue_df[, colSums(!is.na(tissue_df)) > 0]
  
  if (ncol(tissue_df) < 2) next
  
  # Compute correlation matrix
  corr_matrix <- cor(tissue_df, use = "pairwise.complete.obs")
  
  # Save high-resolution PNG
  png(filename = file.path(output_dir, paste0("corrplot_ct_prop_", as.character(t), ".png")),
      width = 8, height = 8, units = "in", res = 300)
  
  # Set very large title font size
  par(cex.main = 3)
  
  # Plot
  corrplot(corr_matrix,
           method = "color",
           type = "full",
           tl.col = "black",
           tl.cex = 1.4,
           tl.srt = 45,
           number.cex = 1.2,
           addCoef.col = "black",
           diag = TRUE,
           title = as.character(t),
           mar = c(0, 0, 5, 0))  # extra top margin
  
  dev.off()
}
