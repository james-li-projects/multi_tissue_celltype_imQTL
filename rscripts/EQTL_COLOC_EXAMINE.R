library(data.table)
library(dplyr)
library(tidyr)

# Define the directory path
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/coloc_output"

# List all files in the directory (modify pattern if needed)
file_list <- list.files(input_dir, full.names = TRUE)

# Read all files with fread and combine into one data.table
coloc_data <- rbindlist(lapply(file_list, fread), fill = TRUE)

# Identify all columns ending in .PP.H4.abf
pp4_cols <- grep("\\.PP\\.H4\\.abf$", names(coloc_data), value = TRUE)

# 1. Create coloc_0.5_all indicator (any PP.H4.abf > 0.5)
coloc_data[, coloc_0.5_all := apply(.SD, 1, function(x) any(x > 0.5, na.rm = TRUE)), .SDcols = pp4_cols]

# 2. For each tissue PP.H4.abf column, create coloc_0.5_{tissue} indicator
for (colname in pp4_cols) {
  tissue <- sub("\\.PP\\.H4\\.abf$", "", colname)
  new_colname <- paste0("coloc_0.5_", tissue)
  coloc_data[, (new_colname) := get(colname) > 0.5]
}

# creating pair_id variable
coloc_data <- coloc_data %>% mutate(pair_id=paste(cpg_id,variant_id,sep="_"))

# separate combination from tissue and celltype
coloc_data <- coloc_data %>% separate(combination,into=c("tissue","celltype"),sep="_",remove=F)

# View result
print(head(coloc_data))

# Examining proportion of coloc overall and for each tissue
coloc_data %>% select(coloc_0.5_all) %>% table()
coloc_data %>% filter(tissue=="colon") %>% select(coloc_0.5_colon) %>% table()
coloc_data %>% filter(tissue=="lung") %>% select(coloc_0.5_lung) %>% table()
coloc_data %>% filter(tissue=="ovary") %>% select(coloc_0.5_ovary) %>% table()
coloc_data %>% filter(tissue=="prostate") %>% select(coloc_0.5_prostate) %>% table()
coloc_data %>% filter(tissue=="wb") %>% select(coloc_0.5_wb) %>% table()

# outputting table
write.table(coloc_data,file="/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/all_eqtl_imqtl_coloc.tsv",quote=F,col.names=T,row.names=F,sep="\t")

#############################################
# CREATING EQTL COLOCALIZATION OUTPUT TABLE #
#############################################
library(dplyr)
library(tidyr)

# Lookup maps
tissue_map <- c("breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
                "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary")
celltype_map <- c("B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell",
                  "NK" = "NK cell", "Neutro" = "Neutrophil", "Mono" = "Monocyte",
                  "IC" = "Immune cells", "DC" = "Dendritic cell", "Epi" = "Epithelial cell",
                  "Endo" = "Endothelial cell", "EndoC" = "Endothelial cell",
                  "MP" = "Macrophage", "Macro" = "Macrophage", "Fib" = "Fibroblast",
                  "Lym" = "Lymphocyte", "Stromal" = "Stromal cell", "Mye" = "Myeloid cell",
                  "EC" = "Endothelial cell", "SM" = "Smooth muscle cell",
                  "LE" = "Luminal epithelium", "BE" = "Basal epithelium")

#########################################
# WRITING OUT SIGNIFICANT COLOC RESULTS #
#########################################
# Apply the mapping
OUTPUT_TABLE_EQTL_COLOC <- coloc_data %>% filter(coloc_0.5_all==TRUE) %>%
  mutate(
    tissue = tissue_map[tissue],
    celltype = celltype_map[celltype]
  ) %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
# outputting table
write.table(OUTPUT_TABLE_EQTL_COLOC,file="/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/OUTPUT_TABLE_EQTL_COLOC.tsv",quote=F,col.names=T,row.names=F,sep="\t")

#################################
# WRITING OUT ALL COLOC RESULTS #
#################################
# Apply the mapping
ALL_OUTPUT_TABLE_EQTL_COLOC <- coloc_data %>%
  mutate(
    tissue = tissue_map[tissue],
    celltype = celltype_map[celltype]
  ) %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
# outputting table
write.table(ALL_OUTPUT_TABLE_EQTL_COLOC,file="/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/ALL_OUTPUT_TABLE_EQTL_COLOC.tsv",quote=F,col.names=T,row.names=F,sep="\t")


###########################################
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(stringr)
library(purrr)

# Step 1: Create summary data
# Define mapping between coloc column names and tissue labels
coloc_cols <- c(
  "coloc_0.5_all" = "All",
  "coloc_0.5_colon" = "Colon",
  "coloc_0.5_lung" = "Lung",
  "coloc_0.5_ovary" = "Ovary",
  "coloc_0.5_prostate" = "Prostate",
  "coloc_0.5_wb" = "Whole Blood"
)

# For each column, filter appropriately and count TRUE/FALSE
pie_data <- imap_dfr(coloc_cols, function(tissue_label, colname) {
  df <- if (colname == "coloc_0.5_all") {
    coloc_data
  } else {
    coloc_data %>% filter(tissue == str_replace(colname, "coloc_0.5_", ""))
  }
  
  df %>%
    count(tissue = tissue_label, coloc_0.5 = .data[[colname]])
}) %>% rename(count=n)
pie_data <- na.omit(pie_data)

# Step 2: Convert to labeled status and compute labels
pie_data <- pie_data %>%
  mutate(
    coloc_status = ifelse(coloc_0.5, "Colocalized", "Not colocalized"),
    coloc_status = factor(coloc_status, levels = c("Colocalized", "Not colocalized"))
  ) %>%
  group_by(tissue) %>%
  mutate(
    percent = count / sum(count) * 100,
    label = paste0(count, " (", format(percent, nsmall = 1, digits = 1), "%)")
  )

# Step 3: Pie chart function with centered slice labels
make_pie <- function(df, tissue_name, show_legend = FALSE) {
  ggplot(df, aes(x = "", y = count, fill = coloc_status)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4.5, color = "black") +
    scale_fill_manual(
      values = c(
        "Colocalized" = "#FA8072",       # salmon
        "Not colocalized" = "#87CEEB"   # sky blue
      ),
      labels = c(
        "Colocalized (PP4 > 0.5)",
        "Not colocalized (PP4 \u2264 0.5)"
      )
    ) +
    labs(title = tissue_name, fill = "Colocalization status") +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16)
    )
}

# Step 4: Generate list of pie charts (only first has legend)
tissues <- unique(pie_data$tissue)
pie_list <- lapply(seq_along(tissues), function(i) {
  t <- tissues[i]
  show_legend <- (i == 1)
  make_pie(pie_data %>% filter(tissue == t), t, show_legend)
})

# Step 5: Combine plots with shared legend and title
final_plot <- wrap_plots(pie_list, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "") &
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# Step 6: Save high-resolution PNG
ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/plots/pie_chart.png",
  plot = final_plot,
  width = 11, height = 9, dpi = 300, units = "in"
)


##################################################
# CREATING HISTOGRAMS OF EACH PP FOR COLOC PAIRS #
##################################################
library(dplyr)
library(tidyr)
library(ggplot2)

tissue_map <- c(
  breast = "Breast", colon = "Colon", lung = "Lung", kidney = "Kidney",
  prostate = "Prostate", wb = "Whole Blood", ovary = "Ovary"
)

pp_colors <- c(
  PP0 = "#1f77b4", PP1 = "#ff7f0e", PP2 = "#2ca02c",
  PP3 = "#d62728", PP4 = "#9467bd"
)

pp_long_df <- OUTPUT_TABLE_EQTL_COLOC %>%
  select(imqtl_id, matches("PP\\.H[0-4]\\.abf$")) %>%
  pivot_longer(
    cols = -imqtl_id,
    names_to = c("tissue", "PP_Hypothesis"),
    names_pattern = "([^.]+)\\.PP\\.(H[0-4])\\.abf",
    values_to = "Posterior_Probability"
  ) %>%
  mutate(
    tissue = tissue_map[tissue] %||% tissue,
    PP_Hypothesis = paste0("PP", substring(PP_Hypothesis, 2))
  ) %>%
  group_by(imqtl_id, tissue) %>%
  filter(any(PP_Hypothesis == "PP4" & Posterior_Probability > 0.5)) %>%
  ungroup()

ggplot(pp_long_df, aes(x = Posterior_Probability, fill = PP_Hypothesis)) +
  geom_histogram(bins = 100, color = "black", alpha = 0.7) +
  scale_fill_manual(values = pp_colors) +
  facet_wrap(~ PP_Hypothesis, ncol = 1, scales = "free_y") +
  labs(x = "Posterior Probability", y = "Count") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> pp_hist_plot

ggsave(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/PP_HIST.png",
  pp_hist_plot, width = 6, height = 8, dpi = 300
)



#######################################
############ PLOT BAR PLOT ############
#######################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Define celltype map
celltype_map <- c(
  "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
  "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
  "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
  "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
  "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
  "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
  "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
)

# Clean up tissue names
coloc_data_label <- coloc_data %>%
  mutate(tissue = recode(tissue,
                         colon = "Colon",
                         lung = "Lung",
                         prostate = "Prostate",
                         ovary = "Ovary",
                         wb = "Whole Blood"))

# Extract celltype from combination and map it
coloc_data_label <- coloc_data_label %>%
  mutate(celltype = str_extract(combination, "[^_]+$")) %>%  # Extract the last segment of combination
  mutate(celltype_display = recode(celltype, !!!celltype_map))  # Map using celltype_map

# creating a dichotomous coloc vs. not coloc variable
coloc_data_label <- coloc_data_label %>% mutate(Category=ifelse(coloc_0.5_all==TRUE,"Colocalizes with an eQTL",ifelse(coloc_0.5_all==FALSE,"No eQTL colocalization",NA)))

# Set tissue order
tissue_order <- rev(c("Prostate", "Ovary", "Lung", "Colon", "Whole Blood"))

# Count per tissue+combo+category
category_counts <- coloc_data_label %>%
  count(tissue, combination, celltype_display, Category) %>%
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) %>%
  mutate(Label = paste0("(", `Colocalizes with an eQTL`, ", ", `No eQTL colocalization`,")"),
         Total = `Colocalizes with an eQTL` + `No eQTL colocalization`)

# Order by tissue then Total
category_counts <- category_counts %>%
  mutate(tissue = factor(tissue, levels = tissue_order)) %>%
  arrange(tissue, desc(Total)) %>%
  mutate(combination = factor(combination, levels = unique(combination)))  # preserve order

# Long format for stacked bars
category_counts_long <- category_counts %>%
  pivot_longer(cols = c(`Colocalizes with an eQTL`, `No eQTL colocalization`), names_to = "Category", values_to = "n")

category_counts_long$Category <- factor(category_counts_long$Category, 
                                        levels = c("Colocalizes with an eQTL", "No eQTL colocalization"))

# Plot
max_x <- max(category_counts$Total, na.rm = TRUE)
max_label_length <- max(nchar(category_counts$Label), na.rm = TRUE)

p <- ggplot(category_counts_long, aes(x = n, y = combination, fill = Category)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(data = category_counts, 
            aes(x = Total + max_x * 0.1 + max_label_length * 0.3, y = combination, label = Label), 
            inherit.aes = FALSE, size = 4, hjust = 0, color = "black") + 
  scale_fill_manual(name = "Directional Consistency",
                    values = c("Colocalizes with an eQTL" = "#00A9A5", "No eQTL colocalization" = "#333333")) +
  scale_y_discrete(labels = category_counts$celltype_display) +  # use mapped labels
  theme_classic() +
  labs(x = "Number of imQTLs colocalizing with eQTLs", y = NULL, fill = "Category") +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    legend.title = element_blank()
  ) +
  expand_limits(x = max_x * 1.5 + max_label_length * 2)

# Save
ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/plots/EQTL_COLOC_STACKED_BARPLOT.png", 
       plot = p, width = 6, height = 6, dpi = 300)

