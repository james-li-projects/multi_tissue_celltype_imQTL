#################################################
############## SINGLE TISSUE PLOTS ##############
#################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(data.table)
library(patchwork)

# importing file
bytissue_enrichment_df <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/bytissue_enrichment_df.tsv")

# Output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/plots/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (tissue_name in c("Colon", "Lung", "Whole Blood")) {
  
  plot_df <- bytissue_enrichment_df %>%
    filter(tissue == tissue_name,
           group != "Chromatin States",
           !is.na(beta_imqtl),
           !is.na(se_imqtl),
           !is.na(beta_mqtl),
           !is.na(se_mqtl)) %>%
    mutate(
      lower_imqtl = beta_imqtl - 1.96 * se_imqtl,
      upper_imqtl = beta_imqtl + 1.96 * se_imqtl,
      lower_mqtl = beta_mqtl - 1.96 * se_mqtl,
      upper_mqtl = beta_mqtl + 1.96 * se_mqtl
    ) %>%
    pivot_longer(
      cols = c(beta_imqtl, beta_mqtl, lower_imqtl, lower_mqtl, upper_imqtl, upper_mqtl),
      names_to = c(".value", "type"),
      names_pattern = "(beta|lower|upper)_(imqtl|mqtl)"
    ) %>%
    mutate(
      type = recode(type, imqtl = "imQTL", mqtl = "mQTL"),
      annotation = as.factor(annotation)
    )
  
  element_plots <- lapply(c("CpG sites", "Variants"), function(el_type) {
    df_el <- plot_df %>%
      filter(element == el_type) %>%
      droplevels() %>%
      group_by(group) %>%
      mutate(annotation = fct_rev(fct_drop(annotation))) %>%
      ungroup()
    
    # Stripe shading
    strip_rects <- df_el %>%
      distinct(group, annotation) %>%
      group_by(group) %>%
      mutate(
        y = rev(row_number()),
        ymin = y - 0.5,
        ymax = y + 0.5,
        shade = rep(c("white", "gray95"), length.out = n())
      ) %>%
      ungroup()
    
    ggplot(df_el, aes(x = beta, y = annotation, color = type)) +
      geom_rect(data = strip_rects, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = shade),
                inherit.aes = FALSE, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(position = position_dodge(width = 0.7)) +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0,
                     position = position_dodge(width = 0.7)) +
      facet_grid(rows = vars(group), scales = "free_y", space = "free_y", switch = "y") +
      scale_color_manual(values = c("imQTL" = "red", "mQTL" = "blue")) +
      scale_fill_identity() +
      labs(
        title = el_type,
        x = "Log Odds Ratio",
        y = NULL,
        color = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
        axis.text.y = element_text(size = 11),
        panel.spacing = unit(1, "lines")
      )
  })
  
  combined_plot <- element_plots[[1]] + element_plots[[2]] + plot_layout(ncol = 2)
  
  ggsave(
    filename = paste0(output_dir, paste0("genomic_context_",tolower(gsub(" ", "_", tissue_name)), ".png")),
    plot = combined_plot,
    width = 10, height = 5, dpi = 300
  )
}


################################################
############## META ANALYSIS PLOT ##############
################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)
library(data.table)

# importing file
meta_enrichment_df <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/meta_enrichment_df.tsv")

# Output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/plots/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare data
plot_df <- meta_enrichment_df %>%
  filter(group != "Chromatin States",
         !is.na(meta_beta_imqtl),
         !is.na(meta_se_imqtl),
         !is.na(meta_beta_mqtl),
         !is.na(meta_se_mqtl)) %>%
  mutate(
    lower_imqtl = meta_beta_imqtl - 1.96 * meta_se_imqtl,
    upper_imqtl = meta_beta_imqtl + 1.96 * meta_se_imqtl,
    lower_mqtl = meta_beta_mqtl - 1.96 * meta_se_mqtl,
    upper_mqtl = meta_beta_mqtl + 1.96 * meta_se_mqtl
  ) %>%
  pivot_longer(
    cols = c(meta_beta_imqtl, meta_beta_mqtl,
             lower_imqtl, lower_mqtl,
             upper_imqtl, upper_mqtl),
    names_to = c(".value", "type"),
    names_pattern = "(meta_beta|lower|upper)_(imqtl|mqtl)"
  ) %>%
  rename(beta = meta_beta) %>%
  mutate(
    type = recode(type, imqtl = "imQTL", mqtl = "mQTL"),
    annotation = as.factor(annotation)
  )

# Generate plot per element
element_plots <- lapply(c("CpG sites", "Variants"), function(el_type) {
  df_el <- plot_df %>%
    filter(element == el_type) %>%
    droplevels() %>%
    group_by(group) %>%
    mutate(annotation = fct_rev(fct_drop(annotation))) %>%
    ungroup()
  
  # Stripe background per group
  strip_rects <- df_el %>%
    distinct(group, annotation) %>%
    group_by(group) %>%
    mutate(
      y = rev(row_number()),
      ymin = y - 0.5,
      ymax = y + 0.5,
      shade = rep(c("white", "gray95"), length.out = n())
    ) %>%
    ungroup()
  
  ggplot(df_el, aes(x = beta, y = annotation, color = type)) +
    geom_rect(data = strip_rects, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = shade),
              inherit.aes = FALSE, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(position = position_dodge(width = 0.7)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0,
                   position = position_dodge(width = 0.7)) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = c("imQTL" = "red", "mQTL" = "blue")) +
    scale_fill_identity() +
    labs(
      title = el_type,
      x = "Log Odds Ratio",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11),
      panel.spacing = unit(1, "lines")
    )
})

# Combine and save
combined_plot <- element_plots[[1]] + element_plots[[2]] + patchwork::plot_layout(ncol = 2)

ggsave(
  filename = paste0(output_dir, "genomic_context_meta.png"),
  plot = combined_plot,
  width = 10, height = 5, dpi = 300
)

