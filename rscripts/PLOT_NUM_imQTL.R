library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#######################################
# extracting numbers of imQTLs
Dataset <- c()
combination <- c()
num_imQTL <- c()

imQTL_feature_df <- data.frame()
imQTL_cpg_df <- data.frame()

for (Dataset in c("GTEx","HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/imQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
  for (combination in file_list) {
    print(combination)
    tmp_df <- fread(combination)
    # filtering for Bonferroni-adjusted p-value<0.05
    tmp_df <- tmp_df %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni"))
    tmp_df <- tmp_df %>% filter(pval_adj_bonf < 0.05)
    num_imQTL <- nrow(tmp_df)
    print(num_imQTL)
    tmp_imQTL_feature_df <- data.frame(Dataset,combination,num_imQTL)
    # assembling a feature DF
    imQTL_feature_df <- rbind(imQTL_feature_df,tmp_imQTL_feature_df)

    # assembling a DF containing all the cpgs
    tmp_imQTL_cpg_df <- tmp_df %>% select(phenotype_id) %>% mutate(combination=combination) %>% mutate(Dataset=Dataset)
    imQTL_cpg_df <- rbind(imQTL_cpg_df,tmp_imQTL_cpg_df)
  }
}

#######################################
# parsing DFs
imQTL_feature_df <- imQTL_feature_df %>% separate(combination,sep="_",into=c("F1","F2","Tissue","celltype_precursor","F3","F4","F5"),remove=F) %>% separate(celltype_precursor, sep= "\\.", into=c("Celltype","F6"))
imQTL_feature_df <- imQTL_feature_df %>% select(Dataset, Tissue, Celltype, num_imQTL, combination) 
imQTL_cpg_df <- imQTL_cpg_df %>% separate(combination,remove=F,sep="_",into=c("F1","F2","Tissue","celltype_precursor","F3","F4","F5")) %>% separate(celltype_precursor, sep= "\\.", into=c("Celltype","F6")) %>% select(Dataset, Tissue, Celltype, phenotype_id, combination) 

#######################################
# making the tissue names in the imQTL_feature_df more beautiful
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="breast"] <- "Breast"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="colon"] <- "Colon"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="lung"] <- "Lung"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="kidney"] <- "Kidney"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="prostate"] <- "Prostate"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="wb"] <- "Whole Blood"
imQTL_feature_df$Tissue[imQTL_feature_df$Tissue=="ovary"] <- "Ovary"
# making the cell type names in the imQTL_feature_df more beautiful
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Basal"] <- "Basal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Luminal"] <- "Luminal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Epi"] <- "Epithelial cell"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="BE"] <- "Basal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="LE"] <- "Luminal epithelium"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="SM"] <- "Smooth muscle cell"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="EC"] <- "Endothelial cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Fat"] <- "Adipocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Fib"] <- "Fibroblast"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Lym"] <- "Lymphocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="MP"] <- "Macrophage"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Macro"] <- "Macrophage"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Mono"] <- "Monocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Mye"] <- "Myeloid cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Gran"] <- "Granulocyte"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Stromal"] <- "Stromal cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Endo"] <- "Endothelial cell"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="B"] <- "B cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="CD4T"] <- "CD4+ T Cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="CD8T"] <- "CD8+ T Cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="NK"] <- "NK cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="Neutro"] <- "Neutrophil"

imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="EndoC"] <- "Endothelial cell"
imQTL_feature_df$Celltype[imQTL_feature_df$Celltype=="IC"] <- "Immune cells"

imQTL_feature_df$imQTL_vec <- paste(imQTL_feature_df$Tissue,imQTL_feature_df$Celltype,sep=" - ")

#######################################
# making the tissue names in the imQTL_cpg_df more beautiful
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="breast"] <- "Breast"
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="colon"] <- "Colon"
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="lung"] <- "Lung"
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="kidney"] <- "Kidney"
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="prostate"] <- "Prostate"
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="wb"] <- "Whole Blood"
imQTL_cpg_df$Tissue[imQTL_cpg_df$Tissue=="ovary"] <- "Ovary"
# making the cell type names in the imQTL_cpg_df more beautiful
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Basal"] <- "Basal epithelium"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Luminal"] <- "Luminal epithelium"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Epi"] <- "Epithelial cell"

imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="BE"] <- "Basal epithelium"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="LE"] <- "Luminal epithelium"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="SM"] <- "Smooth muscle cell"

imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="EC"] <- "Endothelial cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Fat"] <- "Adipocyte"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Fib"] <- "Fibroblast"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Lym"] <- "Lymphocyte"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="MP"] <- "Macrophage"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Macro"] <- "Macrophage"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Mono"] <- "Monocyte"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Mye"] <- "Myeloid cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Gran"] <- "Granulocyte"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Stromal"] <- "Stromal cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Endo"] <- "Endothelial cell"

imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="B"] <- "B cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="CD4T"] <- "CD4+ T Cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="CD8T"] <- "CD8+ T Cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="NK"] <- "NK cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="Neutro"] <- "Neutrophil"

imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="EndoC"] <- "Endothelial cell"
imQTL_cpg_df$Celltype[imQTL_cpg_df$Celltype=="IC"] <- "Immune cells"

imQTL_cpg_df$imQTL_vec <- paste(imQTL_cpg_df$Tissue,imQTL_cpg_df$Celltype,sep=" - ")

#######################################
# filtering out tissue/celltype combinations with no imQTLs and also eosinophils
imQTL_feature_df<-imQTL_feature_df %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Whole Blood")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Breast")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Kidney")) %>% 
  filter(Celltype!="Eosino")
imQTL_cpg_df<-imQTL_cpg_df %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Whole Blood")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Breast")) %>% 
  filter(!(Dataset=="GTEx" & Tissue=="Kidney")) %>% 
  filter(Celltype!="Eosino")

#######################################
# saving DFs
save(imQTL_cpg_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/imQTL_cpg_df.RData")
save(imQTL_feature_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/imQTL_feature_df.RData")

#######################################


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
  "Breast (n=36)" = "turquoise3",
  "Lung (n=190)" = "yellowgreen",
  "Colon (n=189)" = "sienna",
  "Ovary (n=140)" = "pink3",
  "Prostate (n=105)" = "lightgray",
  "Whole Blood (n=1,182)" = "magenta3",
  "Kidney (n=47)" = "turquoise2"
)

####################
# generating plots #
####################
library(ggplot2)
library(scales)

setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis")

png("num_imQTL.png", units = "in", height = 6, width = 4.3, res = 1200)

ggplot(data = modified_imQTL_feature_df %>% filter(num_imQTL > 0), 
       aes(x = imQTL_vec, y = num_imQTL, fill = Tissue)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_y_continuous(trans = 'log10', expand = expansion(mult = c(0.05, 0.15))) +
  coord_flip() +
  labs(x = "Cell Type", y = "Number of mapped imQTLs") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = num_imQTL), hjust = -0.2, size = 3) +
  facet_grid(Dataset ~ ., scales = "free", space = "free") +
  scale_x_discrete(labels = setNames(modified_imQTL_feature_df$Celltype, 
                                     modified_imQTL_feature_df$imQTL_vec)) +
  theme(strip.text.y = element_text(angle = 0),
        panel.border = element_rect(color = "black", fill = NA, size = 0.65),
        legend.position = "bottom",
        plot.margin = margin(5, 20, 5, 5)) +
  scale_fill_manual(values = tissue_colors, guide = guide_legend(nrow = 2))

dev.off()

############################################
# parsing data to create scatterplots of counts of imQTLs
scatterplot_imQTL_feature_df <- modified_imQTL_feature_df %>% mutate(processed_combination_str=gsub("tensorQTL_imQTL_","",gsub(".cis_qtl_top_assoc.txt.gz","",combination)))
scatterplot_imQTL_feature_df$CT_MEAN <- NA
scatterplot_imQTL_feature_df$CT_VAR <- NA
for (j in 1:nrow(scatterplot_imQTL_feature_df)) {
  # extracting tissue and cell type strings
  current_tissue_str <- unlist(str_split(scatterplot_imQTL_feature_df$processed_combination_str[j],pattern="_"))[1]
  current_CT_str <- unlist(str_split(scatterplot_imQTL_feature_df$processed_combination_str[j],pattern="_"))[2]

  # loading cell type proportions
  load(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac/estF.",current_tissue_str,".RData"))
  current_CT_prop <- data.frame(estF.o) %>% dplyr::select(any_of(current_CT_str))
  current_CT_prop_vector <- as.numeric(unlist(current_CT_prop))
  scatterplot_imQTL_feature_df$CT_MEAN[j] <- mean(current_CT_prop_vector)
  scatterplot_imQTL_feature_df$CT_VAR[j] <- var(current_CT_prop_vector)
}
scatterplot_imQTL_feature_df <- scatterplot_imQTL_feature_df %>% rename(`Average proportion`=CT_MEAN) %>% rename(`Variance of proportion`=CT_VAR)
# adding sample size column
scatterplot_imQTL_feature_df <- scatterplot_imQTL_feature_df %>%
  mutate(`Sample Size` = case_when(
    Dataset == "GTEx" & Tissue == "Colon (n=189)" ~ 189,
    Dataset == "GTEx" & Tissue == "Lung (n=190)" ~ 190,
    Dataset == "GTEx" & Tissue == "Ovary (n=140)" ~ 140,
    Dataset == "GTEx" & Tissue == "Prostate (n=105)" ~ 105,
    Dataset == "HEALS" & Tissue == "Whole Blood (n=1,182)" ~ 1182,
    TRUE ~ NA_real_  # Assign NA if no match
  ))





######################################
# create scatterplots with 95% CIs

library(ggplot2)
library(dplyr)

# ==== Helper for p-value formatting ====
format_pval <- function(p) {
  if (p < 0.01) {
    format(p, digits = 2, scientific = TRUE)
  } else {
    format(round(p, 2), nsmall = 2)
  }
}

# ==== Average Proportion Plot ====
plot_df_avg <- scatterplot_imQTL_feature_df %>%
  filter(num_imQTL != 0) %>%
  mutate(
    log_avg = log10(`Average proportion`),
    log_num_imQTL = log10(num_imQTL)
  )

r_val_avg <- round(cor(plot_df_avg$log_avg, plot_df_avg$log_num_imQTL), 2)
model_avg <- lm(log_num_imQTL ~ log_avg, data = plot_df_avg)
p_val_avg_fmt <- format_pval(summary(model_avg)$coefficients[2, "Pr(>|t|)"])

p_avg <- ggplot(plot_df_avg, aes(x = log_avg, y = log_num_imQTL)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
              color = "red", linetype = "dotted", linewidth = 0.6, aes(group = 1)) +
  geom_point(aes(fill = Tissue, size = `Sample Size`), 
             shape = 21, color = "black", alpha = 0.6) +
  annotate("text", x = Inf, y = -Inf,
           label = bquote(italic(r) == .(r_val_avg) ~ "," ~ italic(p) == .(p_val_avg_fmt)),
           hjust = 1.1, vjust = -0.5, size = 7) +
  scale_fill_manual(values = tissue_colors) +
  scale_size_continuous(range = c(6, 18)) +
  labs(
    x = "Mean cell type proportion (log10)",
    y = "Number of mapped imQTLs (log10)"
  ) +
  theme_classic(base_size = 24) +
  theme(legend.position = "bottom",plot.margin = margin(t = 20, r = 30, b = 20, l = 30),legend.spacing.x = unit(0.3, "cm"),legend.text = element_text(size = 11),legend.title = element_text(size = 11)) +
  guides(fill = guide_legend(override.aes = list(size = 6), nrow = 2), size = "none") + 
  ylim(0, 4)

ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/scatterplot_samplesize_numimqtl_avgCT.png",
       plot = p_avg, width = 8, height = 8.25, dpi = 300)


# ==== Variance Proportion Plot ====
plot_df_var <- scatterplot_imQTL_feature_df %>%
  filter(num_imQTL != 0) %>%
  mutate(
    log_var = log10(`Variance of proportion`),
    log_num_imQTL = log10(num_imQTL)
  )

r_val_var <- round(cor(plot_df_var$log_var, plot_df_var$log_num_imQTL), 2)
model_var <- lm(log_num_imQTL ~ log_var, data = plot_df_var)
p_val_var_fmt <- format_pval(summary(model_var)$coefficients[2, "Pr(>|t|)"])

p_var <- ggplot(plot_df_var, aes(x = log_var, y = log_num_imQTL)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
              color = "red", linetype = "dotted", linewidth = 0.6, aes(group = 1)) +
  geom_point(aes(fill = Tissue, size = `Sample Size`), 
             shape = 21, color = "black", alpha = 0.6) +
  annotate("text", x = Inf, y = -Inf,
           label = bquote(italic(r) == .(r_val_var) ~ "," ~ italic(p) == .(p_val_var_fmt)),
           hjust = 1.1, vjust = -0.5, size = 7) +
  scale_fill_manual(values = tissue_colors) +
  scale_size_continuous(range = c(6, 18)) +
  labs(
    x = "Variance of cell type proportion (log10)",
    y = "Number of mapped imQTLs (log10)"
  ) +
  theme_classic(base_size = 24) +
  theme(legend.position = "bottom",plot.margin = margin(t = 20, r = 30, b = 20, l = 30),legend.spacing.x = unit(0.3, "cm"),legend.text = element_text(size = 11),legend.title = element_text(size = 11)) +
  guides(fill = guide_legend(override.aes = list(size = 6), nrow = 2), size = "none") + 
  ylim(0, 4)

ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/scatterplot_samplesize_numimqtl_varCT.png",
       plot = p_var, width = 8, height = 8.25, dpi = 300)

