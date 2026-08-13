library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

################
# MESA INPUTS  #
################

mesa_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/MESA_imQTL_FDR0.05"

mesa_files <- list.files(
  mesa_dir,
  pattern = "table_significant_houseman_.*_imeqtls.txt$",
  full.names = TRUE
)

celltype_map <- c(
  "Bcell" = "wb_B",
  "CD4T"  = "wb_CD4T",
  "CD8T"  = "wb_CD8T",
  "NK"    = "wb_NK",
  "Neu"   = "wb_Neutro"
)

mesa_df <- data.frame()

for (f in mesa_files) {
  
  mesa_celltype <- str_match(
    basename(f),
    "houseman_(.*)_imeqtls"
  )[, 2]
  
  tmp <- fread(f)
  
  tmp <- tmp %>%
    mutate(
      mesa_celltype = mesa_celltype,
      combination = celltype_map[mesa_celltype],
      mesa_variant_id = variant_id,
      variant_id = str_remove(variant_id, "_b38$"),
      variant_id = str_replace_all(variant_id, "_", ":"),
      mesa_b_gi = b_gi,
      mesa_b_gi_se = b_gi_se,
      mesa_pval_gi = pval_gi
    )
  
  mesa_df <- rbind(mesa_df, tmp)
}

mesa_df <- mesa_df %>%
  filter(!is.na(combination))

################
# HEALS INPUTS #
################

pvar <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix.pvar")

afreq <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix_FREQ.afreq") %>%
  filter(ALT_FREQS >= 0.1)

input_dnam_bed <- fread(
  "/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/wb.bed",
  sep = "\t"
)

cov_df <- data.frame(
  t(read.table(
    "/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/processed_covariates_wb.txt",
    sep = "\t"
  ))
)

cov_df <- cov_df %>%
  mutate(Sample_Name = gsub("^X", "", rownames(cov_df)))

write.table(
  data.frame(0, cov_df$Sample_Name),
  file = "/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)

############################
# RUN HEALS VALIDATION     #
############################

all_results <- data.frame()

for (i in seq_len(nrow(mesa_df))) {
  
  target_cpg <- mesa_df$phenotype_id[i]
  lead_variant_raw <- mesa_df$variant_id[i]
  lead_variant_period <- gsub(":", ".", lead_variant_raw)
  tmp_combination <- mesa_df$combination[i]
  
  message("Testing: ", i, " / ", nrow(mesa_df), " | ",
          target_cpg, " | ", lead_variant_raw, " | ", tmp_combination)
  
  #############################
  # Basic CpG / variant checks #
  #############################
  
  tmp_size_cpg <- nrow(
    input_dnam_bed %>%
      filter(phenotype_id == target_cpg) %>%
      head(1)
  )
  
  tmp_size_variant <- nrow(
    afreq %>%
      filter(ID == lead_variant_raw)
  )
  
  if ((tmp_size_cpg + tmp_size_variant) != 2) {
    message("CpG or variant missing; skipping.")
    next
  }
  
  #########################
  # Load HEALS CT vector  #
  #########################
  
  ct_file <- paste0(
    "/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/processed_interactions_",
    tmp_combination,
    ".txt"
  )
  
  ct_df <- read.table(ct_file, sep = "\t") %>%
    rename(Sample_Name = V1, CT = V2)
  
  ########################
  # Export lead genotype #
  ########################
  
  write.table(
    data.frame(lead_variant_raw),
    file = "/gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
  
  system(
    paste0(
      "module load plink/2.0; ",
      "plink2 ",
      "-pfile /gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix ",
      "--extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list ",
      "--keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list ",
      "--maf 0.05 ",
      "--export Av ",
      "--out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_variant"
    )
  )
  
  if (!file.exists("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_variant.traw")) {
    message("PLINK .traw missing; skipping.")
    next
  }
  
  traw <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_variant.traw")
  
  if (nrow(traw) < 1) {
    message("No genotype rows in .traw; skipping.")
    next
  }
  
  traw <- traw %>%
    separate(SNP, into = c("T1", "T2", "T3", "T4"), remove = FALSE) %>%
    mutate(across(11:ncol(.), ~ ifelse(COUNTED != T4, 2 - ., .))) %>%
    select(-c(T1, T2, T3, T4))
  
  variant_df <- data.frame(t(traw))
  newcolnames <- variant_df["SNP", ]
  
  variant_df <- variant_df %>%
    rename_with(~ as.character(newcolnames), everything())
  
  rows_to_remove <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
  
  variant_df <- variant_df %>%
    filter(!rownames(.) %in% rows_to_remove)
  
  variant_df <- variant_df %>%
    mutate(Sample_Name = gsub("^0_", "", rownames(.))) %>%
    select(Sample_Name, everything())
  
  #########################
  # Extract HEALS CpG mval #
  #########################
  
  dnam_bed <- input_dnam_bed %>%
    filter(phenotype_id == target_cpg) %>%
    head(1)
  
  t_dnam_bed <- data.frame(t(dnam_bed))
  
  t_dnam_bed <- t_dnam_bed %>%
    slice(-c(1:4))
  
  cpg_df <- t_dnam_bed %>%
    rename(mval = t.dnam_bed.) %>%
    mutate(Sample_Name = rownames(t_dnam_bed)) %>%
    select(Sample_Name, everything())
  
  #########################
  # Merge regression data #
  #########################
  
  reg_df <- cpg_df %>%
    inner_join(ct_df, by = "Sample_Name") %>%
    inner_join(cov_df, by = "Sample_Name") %>%
    inner_join(variant_df, by = "Sample_Name") %>%
    select(-Sample_Name)
  
  colnames(reg_df) <- gsub(":", ".", colnames(reg_df))
  
  reg_df[] <- lapply(reg_df, function(x) {
    if (is.numeric(x)) {
      x
    } else {
      as.numeric(x)
    }
  })
  
  #########################
  # Define model terms    #
  #########################
  
  chr_indices <- grep("^chr", colnames(reg_df))
  
  lead_variant_index <- which(colnames(reg_df) == lead_variant_period)
  
  if (length(lead_variant_index) != 1) {
    message("Lead variant missing after merge; skipping.")
    next
  }
  
  current_variant_col_name <- colnames(reg_df)[lead_variant_index]
  
  baseline_predictors <- setdiff(
    colnames(reg_df[-c(chr_indices)]),
    "mval"
  )
  
  interaction_term <- paste0(current_variant_col_name, "*CT")
  
  predictors <- c(
    current_variant_col_name,
    baseline_predictors,
    interaction_term
  )
  
  formula <- as.formula(
    paste("mval ~", paste(predictors, collapse = " + "))
  )
  
  #########################
  # Run rank check + lm   #
  #########################
  
  model_matrix <- model.matrix(formula, data = reg_df)
  model_rank <- qr(model_matrix)$rank
  expected_rank <- length(predictors) + 1
  
  if (model_rank < expected_rank) {
    message("Rank deficient model; returning NA.")
    
    heals_beta <- NA
    heals_se <- NA
    heals_t <- NA
    heals_p <- NA
    
  } else {
    
    model <- lm(formula, data = reg_df)
    regression_results <- data.frame(summary(model)$coefficients)
    
    interaction_rowname <- gsub("\\*", ":", interaction_term)
    
    if (!(interaction_rowname %in% rownames(regression_results))) {
      message("Interaction term missing from regression output; skipping.")
      next
    }
    
    heals_beta <- regression_results[interaction_rowname, "Estimate"]
    heals_se <- regression_results[interaction_rowname, "Std..Error"]
    heals_t <- regression_results[interaction_rowname, "t.value"]
    heals_p <- regression_results[interaction_rowname, "Pr...t.."]
  }
  
  #########################
  # Store comparison row  #
  #########################
  
  tmp_result <- data.frame(
    phenotype_id = target_cpg,
    variant_id = lead_variant_raw,
    mesa_variant_id = mesa_df$mesa_variant_id[i],
    mesa_celltype = mesa_df$mesa_celltype[i],
    combination = tmp_combination,
    mesa_b_gi = mesa_df$mesa_b_gi[i],
    mesa_b_gi_se = mesa_df$mesa_b_gi_se[i],
    mesa_pval_gi = mesa_df$mesa_pval_gi[i],
    HEALS_BETA = heals_beta,
    HEALS_SE = heals_se,
    HEALS_T = heals_t,
    HEALS_P = heals_p
  )
  
  all_results <- rbind(all_results, tmp_result)
}

#################
# FINAL OUTPUT  #
#################

all_results <- all_results %>%
  mutate(
    mesa_b_gi = as.numeric(mesa_b_gi),
    mesa_b_gi_se = as.numeric(mesa_b_gi_se),
    mesa_pval_gi = as.numeric(mesa_pval_gi),
    HEALS_BETA = as.numeric(HEALS_BETA),
    HEALS_SE = as.numeric(HEALS_SE),
    HEALS_T = as.numeric(HEALS_T),
    HEALS_P = as.numeric(HEALS_P),
    same_direction = sign(mesa_b_gi) == sign(HEALS_BETA)
  )

out_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/MESA_VALIDATION/MESAFDR0.05_significant_variants_in_HEALS"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save(
  all_results,
  file = file.path(out_dir, "MESA_HEALS_interaction_validation_results.RData")
)


#############################################################
# Validation of MESA FDR0.05 significant variants in HEALS
#############################################################

library(data.table)
library(dplyr)
library(ggplot2)

out_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/MESA_VALIDATION/MESAFDR0.05_significant_variants_in_HEALS/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_map <- data.frame(
  combination = c("wb_B", "wb_CD4T", "wb_CD8T", "wb_NK", "wb_Neutro"),
  mesa_celltype = c("Bcell", "CD4T", "CD8T", "NK", "Neu"),
  plot_name = c("B cell", "CD4 T", "CD8 T", "NK", "Neutrophil")
)

cor_results <- list()

for (i in seq_len(nrow(cell_map))) {
  
  parsed_ct <- cell_map$combination[i]
  mesa_ct   <- cell_map$mesa_celltype[i]
  plot_ct   <- cell_map$plot_name[i]
  
  message("Processing: ", parsed_ct, " / ", mesa_ct)
  
  joined_df <- all_results %>%
    filter(
      combination == parsed_ct,
      mesa_celltype == mesa_ct
    ) %>%
    select(
      phenotype_id, variant_id, mesa_variant_id,
      mesa_celltype, combination,
      mesa_b_gi, mesa_b_gi_se, mesa_pval_gi,
      HEALS_BETA, HEALS_SE, HEALS_T, HEALS_P,
      same_direction
    )
  
  r <- cor(joined_df$HEALS_BETA, joined_df$mesa_b_gi, use = "complete.obs")
  
  message(plot_ct, " correlation: ", round(r, 4), " | n = ", nrow(joined_df))
  
  cor_results[[parsed_ct]] <- data.frame(
    combination = parsed_ct,
    mesa_celltype = mesa_ct,
    n_overlap = nrow(joined_df),
    cor_heals_mesa_b_gi = r
  )
  
  p <- ggplot(joined_df, aes(x = HEALS_BETA, y = mesa_b_gi)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = paste0(plot_ct),
      subtitle = paste0("Pearson r = ", round(r, 3), "; n = ", nrow(joined_df)),
      x = "HEALS interaction effect size",
      y = "MESA interaction effect size (M-values)"
    ) +
    geom_hline(yintercept = 0, color = "gray60", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray60", linewidth = 0.5)
  
  ggsave(
    filename = file.path(out_dir, paste0("MESA_FDR0.05_significant_variants_validated_in_HEALS_", parsed_ct, "_scatter.png")),
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
  file.path(out_dir, "MESA_FDR0.05_significant_variants_validated_in_HEALS_correlations.tsv"),
  sep = "\t"
)