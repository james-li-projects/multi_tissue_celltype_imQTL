#############################################################
# Proportion of HEALS imQTLs colocalizing with any sc-eQTL
# Denominator = unique HEALS imQTLs
#############################################################

library(data.table)
library(dplyr)
library(tidyr)

# Define the directory path
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/custom_prior"

# List coloc output files only
file_list <- list.files(
  input_dir,
  pattern = "^coloc_output_.*\\.tsv$",
  full.names = TRUE
)

# Read all files with fread and combine into one data.table
coloc_data <- rbindlist(lapply(file_list, fread), fill = TRUE)

# Create coloc_0.5_all indicator for each sc-eQTL coloc row
coloc_data[, coloc_0.5_all := PP.H4.abf > 0.5]

# creating pair_id variable
coloc_data <- coloc_data %>%
  mutate(pair_id = paste(cpg_id, imqtl_variant_id, sep = "_"))

# separate parsed_combination from tissue and celltype
coloc_data <- coloc_data %>%
  separate(parsed_combination, into = c("tissue", "celltype"), sep = "_", remove = FALSE)

# View result
print(head(coloc_data))

#############################################################
# Collapse to one row per HEALS imQTL
#############################################################

coloc_data_imqtl <- coloc_data %>%
  group_by(parsed_combination, tissue, celltype, cpg_id, imqtl_variant_id, pair_id) %>%
  summarise(
    coloc_0.5_all = any(coloc_0.5_all, na.rm = TRUE),
    max_PP4 = max(PP.H4.abf, na.rm = TRUE),
    sceqtl_ct_coloc = paste(unique(sceqtl_ct[coloc_0.5_all]), collapse = ";"),
    .groups = "drop"
  )

# Examining proportion of HEALS imQTLs that coloc with any sc-eQTL
coloc_data_imqtl %>% select(coloc_0.5_all) %>% table()

# outputting collapsed table
write.table(
  coloc_data_imqtl,
  file = "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/all_heals_imqtl_sceqtl_coloc_collapsed.tsv",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

#############################################
# CREATING sceQTL COLOCALIZATION OUTPUT TABLE
#############################################

# Lookup maps
tissue_map <- c(
  "wb" = "Whole Blood"
)

celltype_map <- c(
  "B" = "B cell",
  "CD4T" = "CD4+ T Cell",
  "CD8T" = "CD8+ T Cell",
  "NK" = "NK cell",
  "Neutro" = "Neutrophil",
  "Mono" = "Monocyte",
  "Eosino" = "Eosinophil"
)

#########################################
# WRITING OUT SIGNIFICANT COLOC RESULTS
#########################################

OUTPUT_TABLE_SCEQTL_COLOC <- coloc_data_imqtl %>%
  filter(coloc_0.5_all == TRUE) %>%
  mutate(
    tissue = tissue_map[tissue],
    celltype = celltype_map[celltype],
    imqtl_id = paste(parsed_combination, cpg_id, imqtl_variant_id, sep = ".")
  )

write.table(
  OUTPUT_TABLE_SCEQTL_COLOC,
  file = "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/OUTPUT_TABLE_SCEQTL_COLOC.tsv",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

#################################
# WRITING OUT ALL COLOC RESULTS
#################################

ALL_OUTPUT_TABLE_SCEQTL_COLOC <- coloc_data_imqtl %>%
  mutate(
    tissue = tissue_map[tissue],
    celltype = celltype_map[celltype],
    imqtl_id = paste(parsed_combination, cpg_id, imqtl_variant_id, sep = ".")
  )

write.table(
  ALL_OUTPUT_TABLE_SCEQTL_COLOC,
  file = "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/ALL_OUTPUT_TABLE_SCEQTL_COLOC.tsv",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

###########################################
# PIE CHART
###########################################

library(ggplot2)
library(dplyr)
library(tidyr)

# Step 1: Create summary data

pie_data_all <- coloc_data_imqtl %>%
  count(celltype = "All", coloc_0.5 = coloc_0.5_all) %>%
  rename(count = n)

pie_data_celltype <- coloc_data_imqtl %>%
  mutate(celltype_label = celltype_map[celltype]) %>%
  count(celltype = celltype_label, coloc_0.5 = coloc_0.5_all) %>%
  rename(count = n)

pie_data <- bind_rows(pie_data_all, pie_data_celltype)
pie_data <- na.omit(pie_data)

# Step 2: Convert to labeled status and compute labels

pie_data <- pie_data %>%
  mutate(
    coloc_status = ifelse(coloc_0.5, "Colocalized", "Not colocalized"),
    coloc_status = factor(
      coloc_status,
      levels = c("Colocalized", "Not colocalized")
    )
  ) %>%
  group_by(celltype) %>%
  mutate(
    percent = count / sum(count) * 100,
    label = paste0(count, " (", format(percent, nsmall = 1, digits = 1), "%)")
  ) %>%
  ungroup()

# Step 3: Set facet order

celltype_order <- c(
  "All",
  "B cell",
  "CD4+ T Cell",
  "CD8+ T Cell",
  "Monocyte",
  "NK cell",
  "Neutrophil"
)

pie_data$celltype <- factor(
  pie_data$celltype,
  levels = celltype_order[celltype_order %in% unique(pie_data$celltype)]
)

# Step 4: Make faceted pie chart with one shared legend

final_plot <- ggplot(pie_data, aes(x = "", y = count, fill = coloc_status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4.5,
    color = "black"
  ) +
  facet_wrap(~ celltype, ncol = 3, scales = "free") +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Colocalized" = "#FA8072",
      "Not colocalized" = "#87CEEB"
    ),
    labels = c(
      "Colocalized (PP4 > 0.5)",
      "Not colocalized (PP4 \u2264 0.5)"
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

# Step 5: Save high-resolution PNG

ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/custom_prior/pie_chart_heals_imqtls_coloc_with_sceqtls_by_celltype.png",
  plot = final_plot,
  width = 10,
  height = 9,
  dpi = 300,
  units = "in",
  bg = "white"
)