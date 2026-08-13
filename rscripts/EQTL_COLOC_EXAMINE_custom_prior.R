library(data.table)
library(dplyr)
library(tidyr)

# Define the directory path
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/coloc_output_custom_prior"

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
write.table(coloc_data,file="/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/all_eqtl_imqtl_coloc_custom_prior.tsv",quote=F,col.names=T,row.names=F,sep="\t")

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
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/plots/pie_chart_custom_prior.png",
  plot = final_plot,
  width = 11, height = 9, dpi = 300, units = "in"
)
