library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(forcats)
library(patchwork)

# Input files
bytissue_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/bytissue_enrichment_df.tsv"
meta_path     <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/meta_enrichment_df.tsv"

# Output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/plots/"

# Read and prepare data
bytissue_df <- fread(bytissue_path)
meta_df <- fread(meta_path) %>%
  mutate(tissue = "Meta") %>%
  rename(
    beta_imqtl = meta_beta_imqtl, beta_mqtl = meta_beta_mqtl,
    se_imqtl = meta_se_imqtl, se_mqtl = meta_se_mqtl,
    p_imqtl = meta_p_imqtl, p_mqtl = meta_p_mqtl
  )

combined_df <- bind_rows(bytissue_df, meta_df) %>%
  filter(group == "Chromatin States") %>%
  select(tissue, element, annotation, beta_imqtl, p_imqtl, beta_mqtl, p_mqtl) %>%
  pivot_longer(cols = c(beta_imqtl, beta_mqtl, p_imqtl, p_mqtl),
               names_to = c("metric", "type"),
               names_pattern = "(beta|p)_(imqtl|mqtl)") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    stars = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("***", "**", "*", ""), right = FALSE),
    annotation = fct_rev(fct_inorder(annotation)),
    type = recode(type, imqtl = "imQTL", mqtl = "mQTL")
  )

# Shared color scale
beta_range <- combined_df %>%
  summarise(min = min(beta, na.rm = TRUE), max = max(beta, na.rm = TRUE))

# Function to build panel
plot_element_panel <- function(df, element_name, show_legend = TRUE) {
  df_el <- df %>%
    filter(element == element_name) %>%
    mutate(
      is_missing = is.na(beta),
      fill_color = ifelse(is_missing, NA, beta),
      text_label = ifelse(is_missing, "", as.character(stars))
    )
  
  ggplot(df_el, aes(x = type, y = annotation)) +
    geom_tile(aes(fill = fill_color), color = "white") +
    geom_text(aes(label = text_label), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", midpoint = 0,
      limits = c(beta_range$min, beta_range$max),
      na.value = "grey80",
      name = "Log Odds Ratio"
    ) +
    labs(title = element_name, x = NULL, y = NULL) +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "bottom" else "none",
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Generate and save plots per tissue
for (tissue_sel in unique(combined_df$tissue)) {
  df_tissue <- combined_df %>% filter(tissue == tissue_sel)
  
  p_cpg <- plot_element_panel(df_tissue, "CpG sites", show_legend = FALSE)
  p_var <- plot_element_panel(df_tissue, "Variants", show_legend = TRUE)
  
  # Combine horizontally; legend will appear under the "Variants" panel
  final_plot <- p_cpg | p_var
  
  ggsave(
    filename = paste0(output_dir, "chromatin_heatmap_", gsub(" ", "_", tolower(tissue_sel)), ".png"),
    plot = final_plot,
    width = 8, height = 8, dpi = 1200
  )
}
