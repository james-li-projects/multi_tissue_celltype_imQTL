##########################################################
# COMPUTING ENRICHMENT OF EACH ANNOTATION IN EACH TISSUE #
##########################################################
# Initialize result columns
parsed_enrichment_df$beta_imqtl <- NA_real_
parsed_enrichment_df$se_imqtl <- NA_real_
parsed_enrichment_df$p_imqtl <- NA_real_

parsed_enrichment_df$beta_mqtl <- NA_real_
parsed_enrichment_df$se_mqtl <- NA_real_
parsed_enrichment_df$p_mqtl <- NA_real_

for (i in 1:nrow(parsed_enrichment_df)) {
  # ---- imQTL ----
  a <- parsed_enrichment_df$num_intersect_imqtl_elements[i]
  b <- parsed_enrichment_df$num_imqtl_elements[i] - a
  c <- parsed_enrichment_df$num_intersect_background_elements[i]
  d <- parsed_enrichment_df$num_background_elements[i] - c
  
  # Compute logOR, with overflow and Inf/NaN safeguards
  if (b > 0 && c > 0) {
    num <- as.numeric(a) * as.numeric(d)
    den <- as.numeric(b) * as.numeric(c)
    log_or <- log(num / den)
    if (is.finite(log_or)) {
      parsed_enrichment_df$beta_imqtl[i] <- log_or
    }
  }
  
  # SE and P only if no zero counts
  if (all(c(a, b, c, d) > 0)) {
    parsed_enrichment_df$se_imqtl[i] <- sqrt(1/a + 1/b + 1/c + 1/d)
    parsed_enrichment_df$p_imqtl[i] <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
  }
  
  # ---- mQTL ----
  a <- parsed_enrichment_df$num_intersect_mqtl_elements[i]
  b <- parsed_enrichment_df$num_mqtl_elements[i] - a
  c <- parsed_enrichment_df$num_intersect_background_elements[i]
  d <- parsed_enrichment_df$num_background_elements[i] - c
  
  if (b > 0 && c > 0) {
    num <- as.numeric(a) * as.numeric(d)
    den <- as.numeric(b) * as.numeric(c)
    log_or <- log(num / den)
    if (is.finite(log_or)) {
      parsed_enrichment_df$beta_mqtl[i] <- log_or
    }
  }
  
  if (all(c(a, b, c, d) > 0)) {
    parsed_enrichment_df$se_mqtl[i] <- sqrt(1/a + 1/b + 1/c + 1/d)
    parsed_enrichment_df$p_mqtl[i] <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value
  }
  
  if (i %% 10 == 0) cat("Processed row", i, "\n")
}
