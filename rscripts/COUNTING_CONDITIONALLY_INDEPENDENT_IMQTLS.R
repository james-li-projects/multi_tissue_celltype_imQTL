library(data.table)
library(dplyr)
library(stringr)

setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/all_regional_plot_df")

# Set directory containing RData files
rdata_dir <- "./"

# Get list of .RData files
rdata_files <- list.files(rdata_dir, pattern = "\\.RData$", full.names = TRUE)

# Initialize result vector
significant_files <- c()

# Loop through each file
for (file in rdata_files) {
  print(file)
  load(file)  # assumes 'regional_plot_df' is loaded
  
  if (exists("regional_plot_df")) {
    tmp <- subset(regional_plot_df, class == "Interaction Effect (conditional)")
    if (any(tmp$P < 1e-5, na.rm = TRUE)) {
      significant_files <- c(significant_files, basename(file))
    }
    rm(regional_plot_df)
  }
}

conditionally_significant_file_list <- significant_files[!grepl("Eosino",significant_files)]


##############################################
# extracting names and making them beautiful #
##############################################
imQTL_feature_df <- data.frame()
for (combination_index in 1:length(conditionally_significant_file_list)) {
  tmp_combination <- conditionally_significant_file_list[combination_index]
  
  # Extract unprocessed names
  unprocessed_tissue <- str_split(tmp_combination, pattern = "_")[[1]][4]
  unprocessed_celltype <- str_split(tmp_combination, pattern = "_")[[1]][5]
  # Create lookup vectors
  tissue_map <- c(
    "breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
    "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary"
  )
  celltype_map <- c(
    "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
    "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "SM" = "Smooth muscle cell",
    "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
    "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
    "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
    "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell",
    "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
  )
  # Map to beautiful names
  processed_tissue <- tissue_map[[unprocessed_tissue]]
  processed_celltype <- celltype_map[[unprocessed_celltype]]
  processed_tmp_combination <- paste0(processed_tissue," - ",processed_celltype)
  
  # storing results
  tmp_imQTL_feature_df <- data.frame(
    imQTL_vec=processed_tmp_combination,
    Tissue=processed_tissue,
    celltype=processed_celltype
  )
  imQTL_feature_df <- rbind(imQTL_feature_df,tmp_imQTL_feature_df)
}

# get counts
modified_imQTL_feature_df <- imQTL_feature_df %>%
  group_by(imQTL_vec,Tissue,celltype) %>%
  summarise(num_imQTL = n())

################################
# Define the mapping for the Tissue column
tissue_mapping <- c(
  "Lung" = "Lung (n=190)",
  "Colon" = "Colon (n=189)",
  "Ovary" = "Ovary (n=140)",
  "Prostate" = "Prostate (n=105)",
  "Whole Blood" = "Whole Blood (n=1,182)"
)

# Apply the mapping to the Tissue column
modified_imQTL_feature_df <- imQTL_feature_df
modified_imQTL_feature_df$Tissue <- recode(modified_imQTL_feature_df$Tissue, !!!tissue_mapping)

# Define custom colors for specific Tissue values
tissue_colors <- c(
  "Breast" = "turquoise3",
  "Lung" = "yellowgreen",
  "Colon" = "sienna",
  "Ovary" = "pink3",
  "Prostate" = "lightgray",
  "Whole Blood" = "magenta3",
  "Kidney" = "turquoise2"
)

####################
# generating plots #
####################
# barplots for counts of imQTLs
library(ggplot2)
library(scales)
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis")
png("num_imQTL_conditional.png", units = "in", height = 6, width = 6, res = 1200)
ggplot(data = modified_imQTL_feature_df %>% filter(num_imQTL > 0), aes(x = imQTL_vec, y = num_imQTL, fill = Tissue)) +
  geom_bar(stat = "identity", width = 0.75) +
  # scale_y_continuous(trans = 'log10', expand = expansion(mult = c(0.05, 0.15))) + # Add extra space
  coord_flip() +
  labs(x = "Cell Type", y = "Number of mapped imQTLs") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = num_imQTL), hjust = -0.2, size = 2.25) +  # Adjust label position
  theme(strip.text.y = element_text(angle = 0),
        panel.border = element_rect(color = "black", fill = NA, size = 0.65),
        legend.position = "bottom",
        plot.margin = margin(5, 20, 5, 5)) +  # Increase right margin
  scale_fill_manual(values = tissue_colors, guide = guide_legend(nrow = 2)) + 
  guides(fill = guide_legend(nrow = 1))
dev.off()


