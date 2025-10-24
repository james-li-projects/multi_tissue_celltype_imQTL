library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(tibble)
library(rlang)

# ==== USER INPUT ====
imqtl_requests <- CT_specific_imQTL_df %>% mutate(dataset=ifelse(tissue=="wb","HEALS","GTEx"),current_tissue=tissue, target_cpg=phenotype_id, lead_variant_raw=variant_id)

# Load allowed combinations (once)
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")
allowed_combinations <- unique(wide_parsed_imqtl$combination)

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

# Initialize output DF
candidate_ct_specific_imqtl_DF <- data.frame()

# Group by dataset/tissue
grouped_requests <- split(imqtl_requests, paste(imqtl_requests$dataset, imqtl_requests$current_tissue, sep = "_"))

for (group_key in names(grouped_requests)) {
  request_group <- grouped_requests[[group_key]]
  dataset <- unique(request_group$dataset)
  current_tissue <- unique(request_group$current_tissue)
  
  message("Processing dataset = ", dataset, ", tissue = ", current_tissue)
  
  # Load per-tissue data
  pvar <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/", dataset, "/genetic_data/processed_genetic_data_chrprefix.pvar"))
  cov_df <- data.frame(t(read.table(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/", dataset, "/processed_covariates_", current_tissue, ".txt"), sep = "\t")))
  cov_df$Sample_Name <- if (dataset == "GTEx") gsub("\\.", "-", rownames(cov_df)) else gsub("^X", "", rownames(cov_df))
  input_dnam_bed <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/", dataset, "/", current_tissue, ".bed"), sep = "\t")
  
  comb_dir <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/", dataset, "/imQTL/top_assoc")
  tmp_combination_list <- list.files(comb_dir) %>%
    .[grepl(current_tissue, .) & grepl(".cis_qtl_top_assoc.txt.gz", .)] %>%
    gsub("tensorQTL_imQTL_", "", .) %>%
    gsub(".cis_qtl_top_assoc.txt.gz", "", .) %>%
    .[!grepl("Eosino", .)] %>%
    intersect(allowed_combinations)
  
  for (i in seq_len(nrow(request_group))) {
    target_cpg <- request_group$target_cpg[i]
    lead_variant_raw <- request_group$lead_variant_raw[i]
    lead_variant_colon <- lead_variant_raw
    lead_variant_period <- gsub(":", ".", lead_variant_raw)
    alleles <- str_split(lead_variant_raw, ":")[[1]][3:4]
    hom_rec_str <- paste0(alleles[1], "/", alleles[1])
    het_str     <- paste0(alleles[1], "/", alleles[2])
    hom_dom_str <- paste0(alleles[2], "/", alleles[2])
    
    write.table(data.frame(0, cov_df$Sample_Name),
                file = "/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list2",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    
    dnam_bed <- input_dnam_bed %>% filter(phenotype_id == target_cpg) %>% head(1)
    t_dnam_bed <- data.frame(t(dnam_bed))[-c(1:4), , drop = FALSE]
    cpg_df <- t_dnam_bed %>%
      rename(mval = 1) %>%
      mutate(
        Sample_Name = rownames(t_dnam_bed),
        mval = as.numeric(as.character(mval))
      )
    
    # Extract variant ID directly from pvar using exact match
    lead_variant_id <- pvar %>%
      filter(ID == lead_variant_raw) %>%
      pull(ID)
    
    # If not found, skip this request
    if (length(lead_variant_id) == 0) {
      message("Variant not found in pvar: ", lead_variant_raw)
      next
    }
    
    # Write a 1-line variant list file
    write.table(data.frame(ID = lead_variant_id),
                file = "/gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list2",
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    
    system(paste0("module load plink/2.0; plink2 -pfile /gpfs/data/pierce-lab/james.li/imQTL/data/", dataset,
                  "/genetic_data/processed_genetic_data_chrprefix --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list2",
                  " --keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list2 --maf 0.05 --export Av --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window2"))
    
    traw <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window2.traw")
    snp_ids <- gsub(":", ".", traw$SNP)
    genotype_matrix <- as.matrix(traw[, 7:ncol(traw)])
    for (j in seq_len(nrow(traw))) {
      if (traw$COUNTED[j] != traw$ALT[j]) {
        genotype_matrix[j, ] <- 2 - genotype_matrix[j, ]
      }
    }
    variant_df <- as.data.frame(t(genotype_matrix))
    colnames(variant_df) <- snp_ids
    variant_df <- variant_df %>%
      rownames_to_column("Sample_Name") %>%
      mutate(Sample_Name = gsub("^0_", "", Sample_Name))
    
    all_ct_df <- data.frame()
    sig_ct_results <- data.frame()
    
    for (tmp_combination in tmp_combination_list) {
      unprocessed_celltype <- str_split(tmp_combination, "_")[[1]][2]
      if (unprocessed_celltype %in% names(celltype_map)) {
        processed_celltype <- celltype_map[[unprocessed_celltype]]
      } else {
        processed_celltype <- unprocessed_celltype
      }
      processed_tissue <- tissue_map[[current_tissue]] %||% current_tissue
      
      ct_df <- read.table(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/", dataset,
                                 "/processed_interactions_", tmp_combination, ".txt"), sep = "\t") %>%
        rename(Sample_Name = V1, CT = V2)
      
      if (!(lead_variant_period %in% colnames(variant_df))) {
        message("Skipping: variant not in PLINK output - ", lead_variant_period)
        next
      }
      
      reg_df <- cpg_df %>%
        inner_join(ct_df, by = "Sample_Name") %>%
        inner_join(cov_df, by = "Sample_Name") %>%
        inner_join(variant_df, by = "Sample_Name") %>%
        select(Sample_Name, mval, CT, all_of(lead_variant_period)) %>%
        rename(dosage = all_of(lead_variant_period)) %>%
        mutate(
          dosage = round(dosage),
          Genotype = case_when(
            dosage == 0 ~ hom_rec_str,
            dosage == 1 ~ het_str,
            dosage == 2 ~ hom_dom_str,
            TRUE ~ NA_character_
          ),
          `Cell Type Proportion` = ifelse(CT > median(CT),
                                          paste0("Upper 50% (", processed_celltype, ")"),
                                          "All Individuals"),
          Genotype_numeric = as.numeric(as.factor(Genotype))
        ) %>%
        filter(dosage %in% c(0, 1, 2))
      
      # Store for plot
      all_ct_df <- bind_rows(all_ct_df, reg_df)
      
      # Linear model in upper 50%
      upper_df <- reg_df %>% filter(`Cell Type Proportion` == paste0("Upper 50% (", processed_celltype, ")"))
      if (nrow(upper_df) >= 10 && length(unique(upper_df$dosage)) > 1) {
        model <- try(lm(mval ~ dosage, data = upper_df), silent = TRUE)
        if (!inherits(model, "try-error")) {
          pval <- summary(model)$coefficients["dosage", "Pr(>|t|)"]
          sig_ct_results <- bind_rows(sig_ct_results, data.frame(
            celltype = processed_celltype,
            pval = pval
          ))
        }
      }
    }
    
    sig_ct_results <- sig_ct_results %>% filter(pval < 0.05)
    
    if (nrow(sig_ct_results) == 1) {
      candidate_ct_specific_imqtl_DF <- bind_rows(candidate_ct_specific_imqtl_DF, data.frame(
        dataset = dataset,
        tissue = current_tissue,
        lead_variant_raw = lead_variant_raw,
        target_cpg = target_cpg,
        celltype = sig_ct_results$celltype
      ))
      
    } else if (nrow(sig_ct_results) == 2) {
      # Rerun models to extract effect directions
      effect_dirs <- map_dbl(sig_ct_results$celltype, function(ct) {
        upper_df <- all_ct_df %>% filter(`Cell Type Proportion` == paste0("Upper 50% (", ct, ")"))
        if (nrow(upper_df) >= 10 && length(unique(upper_df$dosage)) > 1) {
          model <- try(lm(mval ~ dosage, data = upper_df), silent = TRUE)
          if (!inherits(model, "try-error")) {
            coef(summary(model))["dosage", "Estimate"]
          } else {
            NA
          }
        } else {
          NA
        }
      })
      
      if (all(!is.na(effect_dirs)) && sign(effect_dirs[1]) != sign(effect_dirs[2])) {
        for (k in seq_along(sig_ct_results$celltype)) {
          candidate_ct_specific_imqtl_DF <- bind_rows(candidate_ct_specific_imqtl_DF, data.frame(
            dataset = dataset,
            tissue = current_tissue,
            lead_variant_raw = lead_variant_raw,
            target_cpg = target_cpg,
            celltype = sig_ct_results$celltype[k]
          ))
        }
      }
    }
    
    all_ct_df$`Cell Type Proportion` <- factor(
      all_ct_df$`Cell Type Proportion`,
      levels = c("All Individuals", sort(unique(all_ct_df$`Cell Type Proportion`)[unique(all_ct_df$`Cell Type Proportion`) != "All Individuals"]))
    )
    
    interact_boxplot_p <- ggplot(all_ct_df, aes(x = Genotype, y = mval)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, aes(color = Genotype)) +
      geom_smooth(aes(x = Genotype_numeric, y = mval), method = "lm", se = FALSE, color = "darkgray", linewidth = 1) +
      facet_wrap(~`Cell Type Proportion`, nrow = 1) +
      theme_classic() +
      labs(
        x = "Genotype", y = "DNAm (INT M-value)",
        title = paste(lead_variant_colon, target_cpg, sep = " - ")
      ) +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "none"
      )
    
    out_dir <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/CANDIDATE_CT_SPECIFIC/plots/boxplots/", current_tissue)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(paste0(out_dir, "/", lead_variant_period, "_", target_cpg, ".png"), interact_boxplot_p, width = 20, height = 5)
  }
}

# The resulting data.frame:
candidate_ct_specific_imqtl_DF
write.table(candidate_ct_specific_imqtl_DF,file="/gpfs/data/pierce-lab/james.li/imQTL/output/CANDIDATE_CT_SPECIFIC/tables/candidate_ct_specific_imqtl_DF.tsv",quote=F,row.names=F,col.names=T,sep="\t")
