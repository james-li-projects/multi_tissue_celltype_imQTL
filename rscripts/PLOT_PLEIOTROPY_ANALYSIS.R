library(data.table)

# Define the directory
dir_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/pleiotropy_analysis/tables"

# List all files ending with _variant_counts.tsv
file_list <- list.files(path = dir_path, pattern = "_variant_counts\\.tsv$", full.names = TRUE)

# Read and combine
combined_pleiotropy_df <- rbindlist(lapply(file_list, function(file) {
  # Extract the prefix before _variant_counts.tsv
  combination_name <- sub("_variant_counts\\.tsv$", "", basename(file))
  # Read the file
  df <- fread(file)
  # Add the combination column
  df[, combination := combination_name]
  return(df)
}))

# printing out summary of pleiotropic measures
print(summary(combined_pleiotropy_df$count))
print(sd(combined_pleiotropy_df$count))
print(head(combined_pleiotropy_df%>%arrange(desc(count))))

# Load required library
library(ggplot2)
library(dplyr)

# Define the conversion dictionary
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

# Apply conversion
combined_pleiotropy_df$combination <- conversion_dict[combined_pleiotropy_df$combination]

# Extract tissue from the combination label
combined_pleiotropy_df$tissue <- sapply(strsplit(combined_pleiotropy_df$combination, " - "), `[`, 1)

library(ggplot2)
library(grid)
library(gridExtra)

# Create the plot
p <- ggplot(combined_pleiotropy_df, aes(x = combination, y = count)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_jitter(height = 0, width = 0.2, alpha = 0.3) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of CpGs regulated by lead variants"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 20)  # generous bottom margin
  )

# Convert to grob
g <- ggplotGrob(p)

# Save using grid.draw to avoid clipping
png(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/pleiotropy_analysis/plots/overall_pleiotropy.png",
  width = 12,
  height = 5,
  units = "in",
  res = 300
)
grid::grid.draw(g)
dev.off()
