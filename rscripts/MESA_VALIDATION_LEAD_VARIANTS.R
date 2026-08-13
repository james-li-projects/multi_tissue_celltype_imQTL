#######################################
# importing libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#######################################
# importing whole blood imQTLs
Dataset <- c()
combination <- c()
imqtl_cpg_variant_df <- data.frame()
for (Dataset in c("HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/imQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
  for (combination in file_list) {
    print(combination)
    tmp_df <- fread(combination)
    # filtering for Bonferroni-adjusted p-value<0.05
    tmp_df <- tmp_df %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni"))
    tmp_df <- tmp_df %>% filter(pval_adj_bonf < 0.05)
    
    # assembling a DF containing all the cpgs
    tmp_imqtl_cpg_variant_df <- tmp_df %>% select(phenotype_id,variant_id,b_g,pval_g,b_i,pval_i,b_gi,pval_gi) %>% mutate(combination=combination) %>% mutate(Dataset=Dataset)
    imqtl_cpg_variant_df <- rbind(imqtl_cpg_variant_df,tmp_imqtl_cpg_variant_df)
  }
}
# parsing out chromosome and combination
parsed_imqtl_cpg_variant_df <- imqtl_cpg_variant_df %>% separate(variant_id,remove=F,into=c("chr","pos","a2","a1"),sep="\\:") %>% mutate(chunk_num=as.integer(gsub("chr","",chr))) %>% mutate(parsed_combination=gsub("tensorQTL_imQTL_","",gsub("\\.cis_qtl_top_assoc.txt.gz","",combination))) %>% mutate(input_parquet_filename = paste0(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/",
  Dataset,
  "/imQTL/full_assoc/tensorQTL_imQTL_",
  parsed_combination,
  ".chunk",
  chunk_num,
  ".cis_qtl_pairs.chr",
  chunk_num,
  ".parquet"
)) %>% select(-a2,-a1,-combination,-Dataset)
# filtering out eosinophils
parsed_imqtl_cpg_variant_df <- parsed_imqtl_cpg_variant_df %>% filter(parsed_combination!="wb_Eosino")
# generating a GTEx style variant id
parsed_imqtl_cpg_variant_df <- parsed_imqtl_cpg_variant_df %>% rename(old_variant_id=variant_id) %>% mutate(variant_id = paste0(gsub("\\:","_",old_variant_id),"_b38"))


#############################################################
# Performing validation of interaction effect sizes of lead variants in MESA #
#############################################################
library(data.table)
library(dplyr)
library(ggplot2)

out_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/MESA_VALIDATION/lead_variants/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mesa_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/MESA_imQTL_FDR0.05/"

cell_map <- data.frame(
  parsed_combination = c("wb_B", "wb_CD4T", "wb_CD8T", "wb_NK", "wb_Neutro"),
  mesa_name = c("Bcell", "CD4T", "CD8T", "NK", "Neu"),
  plot_name = c("B cell", "CD4 T", "CD8 T", "NK", "Neutrophil")
)

cor_results <- list()

for (i in seq_len(nrow(cell_map))) {
  
  parsed_ct <- cell_map$parsed_combination[i]
  mesa_ct   <- cell_map$mesa_name[i]
  plot_ct   <- cell_map$plot_name[i]
  
  message("Processing: ", parsed_ct, " / ", mesa_ct)
  
  tmp_imqtl <- parsed_imqtl_cpg_variant_df %>%
    filter(parsed_combination == parsed_ct) %>%
    select(
      phenotype_id, old_variant_id, variant_id,
      b_g, pval_g, b_i, pval_i, b_gi, pval_gi
    ) %>%
    rename(
      heals_b_g = b_g,
      heals_pval_g = pval_g,
      heals_b_i = b_i,
      heals_pval_i = pval_i,
      heals_b_gi = b_gi,
      heals_pval_gi = pval_gi
    )
  
  tmp_mesa <- fread(
    file.path(mesa_dir, paste0("table_significant_houseman_", mesa_ct, "_imeqtls.txt"))
  ) %>%
    select(phenotype_id, variant_id, b_g, pval_g, b_i, pval_i, b_gi, pval_gi) %>%
    rename(
      mesa_b_g = b_g,
      mesa_pval_g = pval_g,
      mesa_b_i = b_i,
      mesa_pval_i = pval_i,
      mesa_b_gi = b_gi,
      mesa_pval_gi = pval_gi
    )
  
  joined_df <- inner_join(
    tmp_imqtl,
    tmp_mesa,
    by = c("phenotype_id", "variant_id")
  )
  
  r <- cor(joined_df$heals_b_gi, joined_df$mesa_b_gi, use = "complete.obs")
  
  message(plot_ct, " correlation: ", round(r, 4), " | n = ", nrow(joined_df))
  
  cor_results[[parsed_ct]] <- data.frame(
    parsed_combination = parsed_ct,
    mesa_name = mesa_ct,
    n_overlap = nrow(joined_df),
    cor_heals_mesa_b_gi = r
  )
  
  p <- ggplot(joined_df, aes(x = heals_b_gi, y = mesa_b_gi)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = paste0(plot_ct),
      subtitle = paste0("Pearson r = ", round(r, 3), "; n = ", nrow(joined_df)),
      x = "HEALS interaction effect size (INT M-values)",
      y = "MESA interaction effect size (M-values)"
    ) + geom_hline(yintercept = 0, color = "gray60", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray60", linewidth = 0.5)
  
  ggsave(
    filename = file.path(out_dir, paste0("MESA_validation_lead_variants_", parsed_ct, "_scatter.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 600
  )
}

cor_results_df <- bind_rows(cor_results)

print(cor_results_df)

fwrite(
  cor_results_df,
  file.path(out_dir, "MESA_validation_lead_variants_correlations.tsv"),
  sep = "\t"
)
