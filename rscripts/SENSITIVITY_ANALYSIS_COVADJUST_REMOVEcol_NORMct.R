library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# ============================================================
# imQTL lead-variant regression runner (NO PLOTTING)
# - For each significant imQTL (Bonferroni<0.05 in top_assoc),
#   fits: mval ~ lead_dosage + CT + lead_dosage:CT + covariates
# - ONLY for the lead variant (no other window variants)
# - Uses pre-built joined_df for ALL covariates + CT fractions
# ============================================================

# ------------------ constants / paths ------------------
MAF_FILTER <- 0.05

tmp_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/tmp"
out_results_path <- "/gpfs/data/pierce-lab/james.li/imQTL/tmp/imQTL_leadvariant_regression_results.tsv.gz"

#dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
#dir.create(dirname(out_results_path), recursive = TRUE, showWarnings = FALSE)

# ------------------ helper: robust numeric conversion ------------------
convert_to_numeric <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) return(col)
    suppressWarnings(as.numeric(col))
  })
  df
}

# ------------------ helper: fit interaction model for lead variant ------------------
# ------------------ helper: fit interaction model for lead variant (NOW also returns main effects) ------------------
fit_lead_imqtl_lm <- function(reg_df, lead_col, covar_cols) {
  fml <- as.formula(
    paste0(
      "mval ~ `", lead_col, "` + CT + `", lead_col, "`:CT",
      if (length(covar_cols) > 0) paste0(" + ", paste(covar_cols, collapse = " + ")) else ""
    )
  )
  
  fit <- tryCatch(stats::lm(fml, data = reg_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  sm <- summary(fit)
  rn <- rownames(sm$coefficients)
  
  # robust term matching
  term_int_candidates <- c(
    paste0("`", lead_col, "`:CT"),
    paste0(lead_col, ":CT"),
    paste0("CT:`", lead_col, "`"),
    paste0("CT:", lead_col)
  )
  term_g_candidates <- c(
    paste0("`", lead_col, "`"),
    lead_col
  )
  term_ct_candidates <- c("CT")
  
  hit_int <- which(rn %in% term_int_candidates)
  hit_g   <- which(rn %in% term_g_candidates)
  hit_ct  <- which(rn %in% term_ct_candidates)
  
  if (length(hit_int) != 1 || length(hit_g) != 1 || length(hit_ct) != 1) return(NULL)
  
  list(
    # genotype main effect
    beta_g = unname(sm$coefficients[hit_g,  "Estimate"]),
    se_g   = unname(sm$coefficients[hit_g,  "Std. Error"]),
    p_g    = unname(sm$coefficients[hit_g,  "Pr(>|t|)"]),
    
    # CT main effect
    beta_ct = unname(sm$coefficients[hit_ct, "Estimate"]),
    se_ct   = unname(sm$coefficients[hit_ct, "Std. Error"]),
    p_ct    = unname(sm$coefficients[hit_ct, "Pr(>|t|)"]),
    
    # interaction effect
    beta_int = unname(sm$coefficients[hit_int, "Estimate"]),
    se_int   = unname(sm$coefficients[hit_int, "Std. Error"]),
    p_int    = unname(sm$coefficients[hit_int, "Pr(>|t|)"]),
    
    n = stats::nobs(fit)
  )
}

# ------------------ main loop ------------------
all_results <- data.table()

gtex_tissue_list <- c("prostate","ovary","lung","colon")
heals_tissue_list <- c("wb")

for (dataset in c("GTEx","HEALS")) {
  
  if (dataset=="GTEx") {
    tissue_list <- gtex_tissue_list
  } else if (dataset=="HEALS") {
    tissue_list <- heals_tissue_list
  }
  
  for (current_tissue in tissue_list) {
    
    # --- importing/parsing covariates/cell types in joined_df ---
    # covariates
    cov_df <- data.frame(t(read.table(
      paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/", dataset, "/processed_covariates_", current_tissue, ".txt"),
      sep = "\t"
    )))
    # cell-type fractions
    load(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/", dataset, "/cell_type_frac/estF.", current_tissue, ".RData"))
    estF.o <- data.frame(estF.o)
    colnames(estF.o) <- paste0("CT_", colnames(estF.o))
    estF.o$Sample_Name <- rownames(estF.o)
    # reformatting ID names based on dataset
    if (dataset=="GTEx") {
      cov_df$Sample_Name <- gsub("\\.", "-", rownames(cov_df))
    } else if (dataset=="HEALS") {
      cov_df$Sample_Name <- gsub("X","",rownames(cov_df))
    } else {
      print("Sample name joining is problematic")
    }
    joined_df <- inner_join(cov_df, estF.o, by = "Sample_Name") %>%
      relocate(any_of("Sample_Name"), .before = everything())
    # drop CT_* columns with >50% zeros
    ct_cols <- grep("^CT_", names(joined_df), value = TRUE)
    ct_remove <- ct_cols[sapply(joined_df[ct_cols], function(x) mean(x == 0, na.rm = TRUE)) > 0.5]
    joined_df <- joined_df %>% select(-all_of(ct_remove))
    # column groups
    dnam_pc_cols <- grep("^PC[0-9]+$", names(joined_df), value = TRUE)   # DNAm PCs only
    gpc_cols     <- grep("^G_PC[0-9]+$", names(joined_df), value = TRUE) # genetic PCs
    ct_cols      <- grep("^CT_", names(joined_df), value = TRUE)
    # correlations (Pearson r + p)
    compute_correlations <- function(df, vars_x, vars_y) {
      expand_grid(x = vars_x, y = vars_y) %>%
        mutate(
          cor     = map2_dbl(x, y, ~ cor(df[[.x]], df[[.y]], use = "pairwise.complete.obs")),
          p_value = map2_dbl(x, y, ~ cor.test(df[[.x]], df[[.y]])$p.value)
        )
    }
    bad_pcs <- union(
      compute_correlations(joined_df, dnam_pc_cols, ct_cols)  %>% filter(p_value < 0.05) %>% pull(x),
      compute_correlations(joined_df, dnam_pc_cols, gpc_cols) %>% filter(p_value < 0.05) %>% pull(x)
    ) %>% sort() %>% unique()
    joined_df <- joined_df %>% select(-all_of(bad_pcs))
   # rank-based INT transforming cell type proportions
    library(RNOmni)
    ct_cols <- grep("^CT_", names(joined_df), value = TRUE)
    print("Rank-based INT transforming CT_ columns using RNOmni::RankNorm")
    setDT(joined_df)
    joined_df[, (ct_cols) := lapply(.SD, RankNorm), .SDcols = ct_cols]
    
    # --- import DNAm bed ---
    input_dnam_bed <- fread(
      paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/", dataset, "/", current_tissue, ".bed"),
      sep = "\t"
    )
    
    # --- sample keep list for plink2 (NOW from joined_df) ---
    keep_path <- file.path(tmp_dir, "tmp_samp.list2")
    write.table(
      data.frame(0, joined_df$Sample_Name),
      file = keep_path,
      quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
    )
    
    # --- find tissue/celltype combinations ---
    setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/", dataset, "/imQTL/top_assoc"))
    
    tmp_combination_list <- list.files()[grepl(current_tissue, list.files()) &
                                           grepl("\\.cis_qtl_top_assoc\\.txt\\.gz$", list.files())]
    tmp_combination_list <- gsub("tensorQTL_imQTL_", "", gsub("\\.cis_qtl_top_assoc\\.txt\\.gz", "", tmp_combination_list))
    tmp_combination_list <- tmp_combination_list[!grepl("Eosino", tmp_combination_list)]  # remove eosinophils
    
    for (tmp_combination in tmp_combination_list) {
      
      message("Dataset: ", dataset, " | Combination: ", tmp_combination)
      
      unprocessed_tissue   <- str_split(tmp_combination, pattern = "_")[[1]][1]
      unprocessed_celltype <- str_split(tmp_combination, pattern = "_")[[1]][2]
      
      tissue_map <- c(
        "breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
        "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary"
      )
      celltype_map <- c(
        "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
        "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "SM" = "Smooth muscle cell",
        "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
        "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
        "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
        "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell",
        "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
      )
      
      processed_tissue   <- tissue_map[[unprocessed_tissue]]
      processed_celltype <- celltype_map[[unprocessed_celltype]]
      
      # --- CT column name in joined_df (e.g., CT_Fib, CT_EndoC, CT_IC) ---
      ct_col <- paste0("CT_", unprocessed_celltype)
      if (!(ct_col %in% colnames(joined_df))) {
        message("  -> Missing ", ct_col, " in joined_df; skipping combination.")
        next
      }
      
      # --- import top assoc results, filter Bonferroni ---
      top_assoc_path <- paste0(
        "/gpfs/data/pierce-lab/james.li/imQTL/output/", dataset,
        "/imQTL/top_assoc/tensorQTL_imQTL_", tmp_combination,
        ".cis_qtl_top_assoc.txt.gz"
      )
      
      tmp_imqtl_df <- fread(top_assoc_path) %>%
        mutate(pval_adj_bonf = p.adjust(pval_emt, method = "bonferroni")) %>%
        filter(pval_adj_bonf < 0.05)
      
      if (nrow(tmp_imqtl_df) < 1) {
        message("  -> No imQTLs at Bonferroni<0.05")
        next
      }
      
      # --- iterate significant imQTLs ---
      for (b in seq_len(nrow(tmp_imqtl_df))) {
        
        target_cpg       <- tmp_imqtl_df$phenotype_id[b]
        lead_variant_raw <- tmp_imqtl_df$variant_id[b]
        lead_variant_period <- gsub(":", ".", lead_variant_raw)
        
        # --- subset DNAm bed to CpG, transpose to sample-wise df ---
        dnam_bed <- input_dnam_bed %>% filter(phenotype_id == target_cpg) %>% head(1)
        if (nrow(dnam_bed) < 1) next
        
        t_dnam_bed <- data.frame(t(dnam_bed)) %>% slice(-c(1:4))
        cpg_df <- t_dnam_bed %>%
          rename(mval = t.dnam_bed.) %>%
          mutate(Sample_Name = rownames(t_dnam_bed)) %>%
          select(Sample_Name, everything())
        
        cpg_df$Sample_Name <- gsub("\\.", "-", cpg_df$Sample_Name)
        
        # --- plink2 extract ONLY the lead variant ---
        extract_path <- file.path(tmp_dir, "lead_variant.list2")
        write.table(
          data.frame(ID = lead_variant_raw),
          file = extract_path,
          quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
        )
        
        plink_prefix <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/", dataset, "/genetic_data/processed_genetic_data_chrprefix")
        out_prefix   <- file.path(tmp_dir, "tmp_lead")
        
        cmd <- paste0(
          "module load plink/2.0; ",
          "plink2 -pfile ", plink_prefix,
          " --extract ", extract_path,
          " --keep ", keep_path,
          " --maf ", MAF_FILTER,
          " --export Av",
          " --out ", out_prefix
        )
        system(cmd)
        
        traw_path <- paste0(out_prefix, ".traw")
        if (!file.exists(traw_path)) next
        
        traw <- fread(traw_path) %>%
          separate(SNP, into = c("T1", "T2", "T3", "T4"), remove = FALSE) %>%
          mutate(across(11:ncol(.), ~ ifelse(COUNTED != T4, 2 - ., .))) %>%
          select(-c(paste0(rep("T", 4), 1:4)))
        
        # --- build sample-wise dosage df (single variant) ---
        variant_df <- data.frame(t(traw))
        newcolnames <- variant_df["SNP", ]
        variant_df <- variant_df %>% rename_with(~ as.character(newcolnames), everything())
        
        rows_to_remove <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
        variant_df <- variant_df %>% filter(!rownames(.) %in% rows_to_remove)
        
        variant_df <- variant_df %>%
          mutate(Sample_Name = gsub("^0_", "", rownames(.))) %>%
          select(Sample_Name, everything())
        
        # --- join into regression df (NOW using joined_df) ---
        # Make a per-combination cov/CT table:
        combo_cov_df <- joined_df %>%
          transmute(
            Sample_Name = Sample_Name,
            CT = .data[[ct_col]],
            across(-all_of(c("Sample_Name", ct_col)))
          )
        
        # Identify other CT columns
        ct_cols <- grep("^CT_", names(combo_cov_df), value = TRUE)
        
        # Drop the last other CT column to avoid collinearity
        combo_cov_df <- combo_cov_df %>%
          select(-all_of(tail(ct_cols, 1))) %>%
          rename_with(~ sub("^CT_", "", .x), starts_with("CT_"))
        
        reg_df <- cpg_df %>%
          inner_join(combo_cov_df, by = "Sample_Name") %>%
          inner_join(variant_df, by = "Sample_Name") %>%
          select(-Sample_Name)
        
        colnames(reg_df) <- gsub(":", ".", colnames(reg_df))
        reg_df <- convert_to_numeric(reg_df)
        
        if (!(lead_variant_period %in% colnames(reg_df))) next
        
        covar_cols <- setdiff(colnames(reg_df), c("mval", "CT", lead_variant_period))
        
        lead_fit <- fit_lead_imqtl_lm(reg_df, lead_variant_period, covar_cols)
        if (is.null(lead_fit)) next
        
        all_results <- rbind(
          all_results,
          data.table(
            dataset = dataset,
            tissue_raw = unprocessed_tissue,
            celltype_raw = unprocessed_celltype,
            tissue = processed_tissue,
            celltype = processed_celltype,
            combination = tmp_combination,
            cpg_id = target_cpg,
            lead_variant_id = lead_variant_raw,
            lead_variant_id_period = lead_variant_period,
            ct_col_used = ct_col,
            
            # NEW: main effects
            beta_g  = lead_fit$beta_g,
            se_g    = lead_fit$se_g,
            p_g     = lead_fit$p_g,
            
            beta_ct = lead_fit$beta_ct,
            se_ct   = lead_fit$se_ct,
            p_ct    = lead_fit$p_ct,
            
            # interaction
            beta_int = lead_fit$beta_int,
            se_int   = lead_fit$se_int,
            p_int    = lead_fit$p_int,
            
            n = lead_fit$n
          ),
          fill = TRUE
        )
      }
    }
  }
}

# ------------------ write results ------------------
if (nrow(all_results) > 0) {
  all_results[, p_adj_bh := p.adjust(p_int, method = "BH")]
  all_results[, p_adj_bonf := p.adjust(p_int, method = "bonferroni")]
}

head(all_results, 1)


write.table(all_results,file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/all_results_REMOVEcol_NORMct.tsv",quote=F,row.names=F,col.names=T,sep="\t")
# fwrite(all_results, out_results_path, sep = "\t")
# message("Wrote: ", out_results_path)