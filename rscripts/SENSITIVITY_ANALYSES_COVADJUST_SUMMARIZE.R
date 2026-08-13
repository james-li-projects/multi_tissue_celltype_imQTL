library(data.table)
library(dplyr)

# importing values from this study's primary cell type imQTL mapping approach
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")

# importing values from sensitivity analyses approaches
files <- c(
  "REMOVEcol_NORMct" = "/gpfs/data/pierce-lab/james.li/imQTL/tmp/all_results_REMOVEcol_NORMct.tsv",
  "ORIGINALcol_RAWct"    = "/gpfs/data/pierce-lab/james.li/imQTL/tmp/all_results_ORIGINALcol_RAWct.tsv"
)

results <- data.table()

for (tag in names(files)) {
  
  cat("Processing:", tag, "\n")
  # importing current new results file
  all_results <- fread(files[tag])
  
  # joining with previous results
  hi1 <- wide_parsed_imqtl %>% select(
    variant_id, phenotype_id, combination,
    Estimate_InteractionEffect, SE_InteractionEffect, P_InteractionEffect,
    Estimate_MainEffect, SE_MainEffect, P_MainEffect
  )
  hi2 <- all_results %>%
    rename(variant_id = lead_variant_id, phenotype_id = cpg_id) %>%
    select(variant_id, phenotype_id, combination, beta_int, se_int, p_int, beta_g, se_g, p_g)
  
  dt <- inner_join(hi1, hi2)
  
  # Force numeric (important)
  dt$Estimate_InteractionEffect <- as.numeric(dt$Estimate_InteractionEffect)
  dt$beta_int <- as.numeric(dt$beta_int)
  dt$Estimate_MainEffect <- as.numeric(dt$Estimate_MainEffect)
  dt$beta_g <- as.numeric(dt$beta_g)
  
  dt$P_InteractionEffect <- as.numeric(dt$P_InteractionEffect)
  dt$p_int <- as.numeric(dt$p_int)
  dt$P_MainEffect <- as.numeric(dt$P_MainEffect)
  dt$p_g <- as.numeric(dt$p_g)
  
  # Remove p-values <= 0 (avoid log10(0))
  dt$P_InteractionEffect[dt$P_InteractionEffect <= 0] <- NA
  dt$p_int[dt$p_int <= 0] <- NA
  dt$P_MainEffect[dt$P_MainEffect <= 0] <- NA
  dt$p_g[dt$p_g <= 0] <- NA
  
  # ---- Correlations ----
  cor_beta_inter <- cor(dt$Estimate_InteractionEffect, dt$beta_int, use = "pairwise.complete.obs")
  cor_beta_main  <- cor(dt$Estimate_MainEffect, dt$beta_g, use = "pairwise.complete.obs")
  
  cor_logp_inter <- cor(-log10(dt$P_InteractionEffect), -log10(dt$p_int), use = "pairwise.complete.obs")
  cor_logp_main  <- cor(-log10(dt$P_MainEffect), -log10(dt$p_g), use = "pairwise.complete.obs")
  
  # ---- p-value summaries (min, Q1, median, mean, Q3, max) ----
  inter_p_summ <- quantile(dt$p_int, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  main_p_summ  <- quantile(dt$p_g,   probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  
  inter_p_mean <- mean(dt$p_int, na.rm = TRUE)
  main_p_mean  <- mean(dt$p_g,   na.rm = TRUE)
  
  results <- rbind(
    results,
    data.table(
      file = tag,
      cor_beta_inter = cor_beta_inter,
      cor_beta_main = cor_beta_main,
      cor_logp_inter = cor_logp_inter,
      cor_logp_main = cor_logp_main,
      
      # interaction term p-value distribution (from sensitivity file: p_int)
      p_int_min  = inter_p_summ[1],
      p_int_q1   = inter_p_summ[2],
      p_int_med  = inter_p_summ[3],
      p_int_mean = inter_p_mean,
      p_int_q3   = inter_p_summ[4],
      p_int_max  = inter_p_summ[5],
      
      # main effect p-value distribution (from sensitivity file: p_g)
      p_g_min  = main_p_summ[1],
      p_g_q1   = main_p_summ[2],
      p_g_med  = main_p_summ[3],
      p_g_mean = main_p_mean,
      p_g_q3   = main_p_summ[4],
      p_g_max  = main_p_summ[5]
    )
  )
}

# ---- ADD: p-value summaries from PRIMARY imQTL mapping ----
primary_p_int <- as.numeric(wide_parsed_imqtl$P_InteractionEffect)
primary_p_g   <- as.numeric(wide_parsed_imqtl$P_MainEffect)

primary_p_int[primary_p_int <= 0] <- NA
primary_p_g[primary_p_g <= 0] <- NA

primary_inter <- quantile(primary_p_int, probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE, names = FALSE)
primary_main  <- quantile(primary_p_g,   probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE, names = FALSE)

primary_inter_mean <- mean(primary_p_int, na.rm = TRUE)
primary_main_mean  <- mean(primary_p_g,   na.rm = TRUE)

results <- rbind(
  results,
  data.table(
    file = "PRIMARY_imQTL",
    
    cor_beta_inter = NA,
    cor_beta_main = NA,
    cor_logp_inter = NA,
    cor_logp_main = NA,
    
    p_int_min  = primary_inter[1],
    p_int_q1   = primary_inter[2],
    p_int_med  = primary_inter[3],
    p_int_mean = primary_inter_mean,
    p_int_q3   = primary_inter[4],
    p_int_max  = primary_inter[5],
    
    p_g_min  = primary_main[1],
    p_g_q1   = primary_main[2],
    p_g_med  = primary_main[3],
    p_g_mean = primary_main_mean,
    p_g_q3   = primary_main[4],
    p_g_max  = primary_main[5]
  )
)

print(results)

write.table(
  x = data.frame(t(results)) %>% select(X3, X1, X2),
  file = "/gpfs/data/pierce-lab/james.li/imQTL/output/SENSITIVITY_ANALYSES/sensitivity_analyses.tsv",
  quote = FALSE,
  row.names = TRUE,
  col.names = FALSE,
  sep = "\t"
)