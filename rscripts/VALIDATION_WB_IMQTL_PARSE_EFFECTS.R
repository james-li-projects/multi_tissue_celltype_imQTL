########################################
# plotting scatter plots of all imQTLs #
########################################
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)

# Load data
load("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/all_imqtl_results.RData")

# Rename combinations
all_imqtl_results <- all_imqtl_results %>%
  mutate(combination = recode(combination,
                              "wb_B" = "B cell",
                              "wb_CD4T" = "CD4+ T cell",
                              "wb_CD8T" = "CD8+ T cell",
                              "wb_Mono" = "Monocyte",
                              "wb_NK" = "NK cell",
                              "wb_Neutro" = "Neutrophil"
  ))

# creating a paired DF
### for each cell type
all_imqtl_results <- all_imqtl_results %>% select(-SE,-T,-class)
all_imqtl_results_GTEx <- all_imqtl_results %>% filter(dataset=="GTEx") %>% select(-dataset) %>% rename(GTEx=BETA,P_GTEx=P)
all_imqtl_results_HEALS <- all_imqtl_results %>% filter(dataset=="HEALS") %>% select(-dataset) %>% rename(HEALS=BETA,P_HEALS=P)
paired_imqtl_df_part1 <- inner_join(all_imqtl_results_HEALS,all_imqtl_results_GTEx,by=c("variant_id","phenotype_id","combination")) %>% na.omit()
### for all cell types
paired_imqtl_df_part2 <- paired_imqtl_df_part1 %>% mutate(combination="All cell types")
### combined
paired_imqtl_df <- rbind(paired_imqtl_df_part1,paired_imqtl_df_part2)

# Factor for facet ordering
celltype_levels <- c(
  "All cell types", "B cell", "CD4+ T cell", "CD8+ T cell",
  "Monocyte", "NK cell", "Neutrophil"
)
paired_imqtl_df$combination <- factor(paired_imqtl_df$combination, levels = celltype_levels)

# Generate regression labels
regression_labels <- paired_imqtl_df %>%
  group_by(combination) %>%
  do({
    tryCatch({
      model <- lm(HEALS ~ GTEx, data = .)
      summary_model <- summary(model)
      r2 <- summary_model$r.squared
      pval <- coef(summary_model)[2, 4]
      
      r2_label <- sprintf("r2 = %.2f", r2)
      pval_label <- if (!is.na(pval) && pval < 0.05) {
        paste0("p = ", formatC(pval, format = "e", digits = 1))
      } else {
        paste0("p = ", sprintf("%.2f", pval))
      }
      
      tibble(combination = unique(.$combination),
             label = paste(r2_label, pval_label, sep = "\n"))
    }, error = function(e) {
      tibble(combination = unique(.$combination), label = "Insufficient data")
    })
  }) %>%
  ungroup()

# Match factor levels in labels
regression_labels$combination <- factor(regression_labels$combination, levels = celltype_levels)

# Restrict smoothing to panels with enough data
df_for_smooth <- paired_imqtl_df %>%
  group_by(combination) %>%
  filter(n() >= 3, n_distinct(GTEx) > 1, n_distinct(HEALS) > 1) %>%
  ungroup()

# Plot
p <- ggplot(paired_imqtl_df, aes(x = HEALS, y = GTEx)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(data = df_for_smooth, method = "lm", se = TRUE, color = "red") +
  facet_wrap(~ combination, scales = "free", ncol = 3) +
  geom_text(data = regression_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.3, inherit.aes = FALSE, size = 3.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  theme_classic() +
  xlab("HEALS Interaction Estimate") +
  ylab("GTEx Interaction Estimate") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

# Save to file
ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/scatterplot_celltype_imqtl_all.png",
       plot = p, width = 10, height = 6, dpi = 600)



##########################################################
# plotting scatter plots of imQTLs with a p<0.05 in GTEx #
##########################################################
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)

# Load data
load("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/all_imqtl_results.RData")

# Rename combinations
all_imqtl_results <- all_imqtl_results %>%
  mutate(combination = recode(combination,
                              "wb_B" = "B cell",
                              "wb_CD4T" = "CD4+ T cell",
                              "wb_CD8T" = "CD8+ T cell",
                              "wb_Mono" = "Monocyte",
                              "wb_NK" = "NK cell",
                              "wb_Neutro" = "Neutrophil"
  ))

# creating a paired DF
### for each cell type
all_imqtl_results <- all_imqtl_results %>% select(-SE,-T,-class)
all_imqtl_results_GTEx <- all_imqtl_results %>% filter(dataset=="GTEx") %>% select(-dataset) %>% rename(GTEx=BETA,P_GTEx=P)
all_imqtl_results_HEALS <- all_imqtl_results %>% filter(dataset=="HEALS") %>% select(-dataset) %>% rename(HEALS=BETA,P_HEALS=P)
paired_imqtl_df_part1 <- inner_join(all_imqtl_results_HEALS,all_imqtl_results_GTEx,by=c("variant_id","phenotype_id","combination")) %>% na.omit()
### for all cell types
paired_imqtl_df_part2 <- paired_imqtl_df_part1 %>% mutate(combination="All cell types")
### combined
paired_imqtl_df <- rbind(paired_imqtl_df_part1,paired_imqtl_df_part2)

# filtering for GTEx variants with a p-value < 0.05
paired_imqtl_df <- paired_imqtl_df %>% filter(P_GTEx<0.05)

# Factor for facet ordering
celltype_levels <- c(
  "All cell types", "B cell", "CD4+ T cell", "CD8+ T cell",
  "Monocyte", "NK cell", "Neutrophil"
)
paired_imqtl_df$combination <- factor(paired_imqtl_df$combination, levels = celltype_levels)

# Generate regression labels
regression_labels <- paired_imqtl_df %>%
  group_by(combination) %>%
  do({
    tryCatch({
      model <- lm(HEALS ~ GTEx, data = .)
      summary_model <- summary(model)
      r2 <- summary_model$r.squared
      pval <- coef(summary_model)[2, 4]
      
      r2_label <- sprintf("r2 = %.2f", r2)
      pval_label <- if (!is.na(pval) && pval < 0.05) {
        paste0("p = ", formatC(pval, format = "e", digits = 1))
      } else {
        paste0("p = ", sprintf("%.2f", pval))
      }
      
      tibble(combination = unique(.$combination),
             label = paste(r2_label, pval_label, sep = "\n"))
    }, error = function(e) {
      tibble(combination = unique(.$combination), label = "Insufficient data")
    })
  }) %>%
  ungroup()

# Match factor levels in labels
regression_labels$combination <- factor(regression_labels$combination, levels = celltype_levels)

# Restrict smoothing to panels with enough data
df_for_smooth <- paired_imqtl_df %>%
  group_by(combination) %>%
  filter(n() >= 3, n_distinct(GTEx) > 1, n_distinct(HEALS) > 1) %>%
  ungroup()

# Plot
p <- ggplot(paired_imqtl_df, aes(x = HEALS, y = GTEx)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(data = df_for_smooth, method = "lm", se = TRUE, color = "red") +
  facet_wrap(~ combination, scales = "free", ncol = 3) +
  geom_text(data = regression_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.3, inherit.aes = FALSE, size = 3.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  theme_classic() +
  xlab("HEALS Interaction Estimate") +
  ylab("GTEx Interaction Estimate") +
  theme(
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

# Save to file
ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/scatterplot_celltype_imqtl_p0.05.png",
       plot = p, width = 10, height = 6, dpi = 600)


##################################
# Creating validation pie charts #
##################################
library(dplyr)
library(ggplot2)
library(broom)
library(tidyr)
library(data.table)
library(tidyr)
library(stringr)

# Load data
load("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/all_imqtl_results.RData")

# Rename combinations
all_imqtl_results <- all_imqtl_results %>%
  mutate(combination = recode(combination,
                              "wb_B" = "B cell",
                              "wb_CD4T" = "CD4+ T cell",
                              "wb_CD8T" = "CD8+ T cell",
                              "wb_Mono" = "Monocyte",
                              "wb_NK" = "NK cell",
                              "wb_Neutro" = "Neutrophil"
  ))

# creating a paired DF
### for each cell type
all_imqtl_results <- all_imqtl_results %>% select(-SE,-T,-class)
all_imqtl_results_GTEx <- all_imqtl_results %>% filter(dataset=="GTEx") %>% select(-dataset) %>% rename(b_gi_GTEx=BETA,P_GTEx=P)
all_imqtl_results_HEALS <- all_imqtl_results %>% filter(dataset=="HEALS") %>% select(-dataset) %>% rename(b_gi_HEALS=BETA,P_HEALS=P)
paired_imqtl_df_part1 <- inner_join(all_imqtl_results_HEALS,all_imqtl_results_GTEx,by=c("variant_id","phenotype_id","combination")) %>% na.omit()
### for all cell types
paired_imqtl_df_part2 <- paired_imqtl_df_part1 %>% mutate(combination="All cell types")
### combined
paired_imqtl_df <- rbind(paired_imqtl_df_part1,paired_imqtl_df_part2)

# further parsing this combined DF
paired_imqtl_df <- paired_imqtl_df %>% mutate(variant_id=gsub("\\.",":",variant_id))
paired_imqtl_df <- paired_imqtl_df %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_")) 
# adding validation status
### validation for just sign
paired_imqtl_df_only_sign_validation <- paired_imqtl_df %>%
  mutate(validated_status = sign(b_gi_GTEx) == sign(b_gi_HEALS))
paired_imqtl_df_only_sign_validation <- paired_imqtl_df_only_sign_validation %>% rename(celltype=combination) %>% select(pair_id,b_gi_HEALS,b_gi_GTEx,validated_status,celltype)
### official validation (sign + p<0.05)
paired_imqtl_df <- paired_imqtl_df %>%
  mutate(validated_status = P_GTEx < 0.05 & sign(b_gi_GTEx) == sign(b_gi_HEALS))
paired_imqtl_df <- paired_imqtl_df %>% rename(celltype=combination) %>% select(pair_id,b_gi_HEALS,b_gi_GTEx,validated_status,celltype)


# extract distances for these whole blood HEALS imQTLs
wb_distance_imQTL_cpg_variant_df <- data.frame()
for (Dataset in c("HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/imQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
  file_list <- file_list[grepl("wb",file_list)]
  file_list <- file_list[!grepl("Eosino",file_list)]
  
  for (combination in file_list) {
    print(combination)
    tmp_df <- fread(combination)
    # filtering for Bonferroni-adjusted p-value<0.05
    tmp_df <- tmp_df %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni"))
    tmp_df <- tmp_df %>% filter(pval_adj_bonf < 0.05)
    
    # assembling a DF containing all the cpgs
    tmp_wb_distance_imQTL_cpg_variant_df <- tmp_df %>% select(phenotype_id,variant_id,start_distance,end_distance) %>% mutate(combination=combination) %>% mutate(Dataset=Dataset)
    wb_distance_imQTL_cpg_variant_df <- rbind(wb_distance_imQTL_cpg_variant_df,tmp_wb_distance_imQTL_cpg_variant_df)
  }
}
# Define mapping from file names to cell type names
combination_to_celltype <- c(
  "tensorQTL_imQTL_wb_B.cis_qtl_top_assoc.txt.gz" = "B cell",
  "tensorQTL_imQTL_wb_CD4T.cis_qtl_top_assoc.txt.gz" = "CD4+ T cell",
  "tensorQTL_imQTL_wb_CD8T.cis_qtl_top_assoc.txt.gz" = "CD8+ T cell",
  "tensorQTL_imQTL_wb_Mono.cis_qtl_top_assoc.txt.gz" = "Monocyte",
  "tensorQTL_imQTL_wb_NK.cis_qtl_top_assoc.txt.gz" = "NK cell",
  "tensorQTL_imQTL_wb_Neutro.cis_qtl_top_assoc.txt.gz" = "Neutrophil"
)
# Apply mapping
wb_distance_imQTL_cpg_variant_df <- wb_distance_imQTL_cpg_variant_df %>%
  mutate(celltype = combination_to_celltype[combination]) 

# parsing wb_distance_imQTL_cpg_variant_df
wb_distance_imQTL_cpg_variant_df <- wb_distance_imQTL_cpg_variant_df %>% mutate(pair_id = paste(phenotype_id, variant_id,sep="_"),mean_distance=(start_distance+end_distance)/2) %>% mutate(abs_distance=abs(mean_distance))
wb_distance_imQTL_cpg_variant_df <- wb_distance_imQTL_cpg_variant_df %>% select(start_distance,end_distance,pair_id,mean_distance,abs_distance,celltype)

################################
# joining validation and distance DFs
all_joined_df <- inner_join(wb_distance_imQTL_cpg_variant_df,paired_imqtl_df,by=c("pair_id", "celltype")) %>% select(start_distance,end_distance,pair_id,b_gi_HEALS,b_gi_GTEx,validated_status,mean_distance,abs_distance,celltype)

# also joining the validation based solely on sign and distance DFs
all_joined_df_only_sign_validation <- inner_join(wb_distance_imQTL_cpg_variant_df,paired_imqtl_df_only_sign_validation,by=c("pair_id", "celltype")) %>% select(start_distance,end_distance,pair_id,b_gi_HEALS,b_gi_GTEx,validated_status,mean_distance,abs_distance,celltype)
################################





############################
# VALIDATION DISTANCE PLOT #
############################
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

# just for reference, examinine distance by validation status
all_joined_df %>%
  group_by(validated_status) %>%
  summarise(
    count = n(),
    min = min(abs_distance, na.rm = TRUE),
    q1 = quantile(abs_distance, 0.25, na.rm = TRUE),
    median = median(abs_distance, na.rm = TRUE),
    mean = mean(abs_distance, na.rm = TRUE),
    q3 = quantile(abs_distance, 0.75, na.rm = TRUE),
    max = max(abs_distance, na.rm = TRUE),
    sd = sd(abs_distance, na.rm = TRUE)
  )

# Colorblind-friendly palette (reversed)
cb_palette <- c("Not validated" = "#0072B2", "Validated" = "#E69F00")

# Prepare data
all_joined_df_binned <- all_joined_df %>%
  filter(mean_distance < 100000 & mean_distance > -100000) %>%
  mutate(
    distance_bin = cut(mean_distance,
                       breaks = seq(floor(min(mean_distance, na.rm = TRUE) / 5000) * 5000,
                                    ceiling(max(mean_distance, na.rm = TRUE) / 5000) * 5000,
                                    by = 5000),
                       right = FALSE, include.lowest = TRUE)
  )

# Count entries per bin and validated status
distance_counts <- all_joined_df_binned %>%
  group_by(distance_bin, validated_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    bin_range = gsub("\\[|\\)|\\]", "", distance_bin)
  ) %>%
  separate(bin_range, into = c("start", "end"), sep = ",", convert = TRUE) %>%
  mutate(
    midpoint_kb = -1*((start + end) / 2 / 1000),  # signed midpoint in kb
    validated_label = ifelse(validated_status, "Validated", "Not validated")
  )

# Plot
p <- ggplot(distance_counts, aes(x = midpoint_kb, y = n, fill = validated_label)) +
  geom_bar(stat = "identity", position = "stack", width = 5, color = "black") +
  scale_fill_manual(values = cb_palette) +
  scale_x_continuous(
    breaks = seq(-100, 100, by = 20)
  ) +
  labs(
    x = "CpG to variant distance (kb)",
    y = "Number of imQTLs",
    fill = ""
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Save high-resolution PNG
ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/distance_validate_hist.png",
       plot = p, width = 6.5, height = 6, dpi = 1200)



#################################################
# VALIDATION PIE CHART PLOT OFFICIAL VALIDATION #
#################################################
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

all_joined_df$celltype <- as.character(all_joined_df$celltype)
all_joined_df$celltype[is.na(all_joined_df$celltype)] <- "Unknown"

preferred_levels <- c("All whole blood imQTLs", "B cell", "CD4+ T cell", "CD8+ T cell", "NK cell", "Neutrophil")

pie_data <- all_joined_df %>%
  group_by(celltype, validated_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(
    percent = n / sum(n),
    validated_label = ifelse(validated_status, "Validated", "Not validated"),
    label = paste0(n, "\n(", formatC(percent * 100, format = "f", digits = 1), "%)")
  )

pie_data_all <- all_joined_df %>%
  mutate(celltype = "All whole blood imQTLs") %>%
  group_by(celltype, validated_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(
    percent = n / sum(n),
    validated_label = ifelse(validated_status, "Validated", "Not validated"),
    label = paste0(n, "\n(", formatC(percent * 100, format = "f", digits = 1), "%)")
  )



# joining all pie chart DFs and removing monocytes since there was too low of a sample size
all_joined_df <- all_joined_df 
pie_data_combined <- bind_rows(pie_data, pie_data_all) %>% filter(celltype!="Monocyte") %>%  mutate(
    celltype = factor(celltype, levels = preferred_levels),
    validated_label = factor(validated_label, levels = c("Validated", "Not validated"))
  )

cb_palette <- c("Validated" = "#E69F00", "Not validated" = "#0072B2")

p <- ggplot(pie_data_combined, aes(x = "", y = percent, fill = validated_label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) +
  facet_wrap(~celltype, nrow = 2) +
  scale_fill_manual(values = cb_palette) +
  labs(fill = "Validation status") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    panel.border = element_blank()
  )

ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/validation_piecharts.png",
       plot = p, width = 12, height = 7, dpi = 600)








##################################################
# VALIDATION PIE CHART PLOT ONLY SIGN VALIDATION #
##################################################
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

all_joined_df_only_sign_validation$celltype <- as.character(all_joined_df_only_sign_validation$celltype)
all_joined_df_only_sign_validation$celltype[is.na(all_joined_df_only_sign_validation$celltype)] <- "Unknown"

preferred_levels <- c("All whole blood imQTLs", "B cell", "CD4+ T cell", "CD8+ T cell", "NK cell", "Neutrophil")

pie_data <- all_joined_df_only_sign_validation %>%
  group_by(celltype, validated_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(
    percent = n / sum(n),
    validated_label = ifelse(validated_status, "Same sign", "Different sign"),
    label = paste0(n, "\n(", formatC(percent * 100, format = "f", digits = 1), "%)")
  )

pie_data_all <- all_joined_df_only_sign_validation %>%
  mutate(celltype = "All whole blood imQTLs") %>%
  group_by(celltype, validated_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(
    percent = n / sum(n),
    validated_label = ifelse(validated_status, "Same sign", "Different sign"),
    label = paste0(n, "\n(", formatC(percent * 100, format = "f", digits = 1), "%)")
  )



# joining all pie chart DFs and removing monocytes since there was too low of a sample size
all_joined_df_only_sign_validation <- all_joined_df_only_sign_validation 
pie_data_combined <- bind_rows(pie_data, pie_data_all) %>% filter(celltype!="Monocyte") %>%  mutate(
  celltype = factor(celltype, levels = preferred_levels),
  validated_label = factor(validated_label, levels = c("Same sign", "Different sign"))
)

cb_palette <- c("Same sign" = "#E69F00", "Different sign" = "#0072B2")

p <- ggplot(pie_data_combined, aes(x = "", y = percent, fill = validated_label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) +
  facet_wrap(~celltype, nrow = 2) +
  scale_fill_manual(values = cb_palette) +
  labs(fill = "Directional consistency in sign of interaction terms") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    panel.border = element_blank()
  )

ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/validation_piecharts_ONLY_SIGN_VALIDATION.png",
       plot = p, width = 12, height = 7, dpi = 600)


