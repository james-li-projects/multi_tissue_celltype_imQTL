library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#############################################################
## SETTING UP LISTS OF ELEMENTS FOR IMQTLS/MQTLS/BACKGROUND #
#############################################################
#######################################
# extracting cpg and variant IDs for imQTLs
imQTL_cpg_variant_df <- data.frame()
for (Dataset in c("GTEx","HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/imQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
  for (combination in file_list) {
    print(combination)
    tmp_df <- fread(combination)
    
    # filtering for Bonferroni-adjusted p-value<0.05
    tmp_df <- tmp_df %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni"))
    tmp_df <- tmp_df %>% filter(pval_adj_bonf < 0.05)
    
    # assembling a DF containing all the cpg and variant IDs
    tmp_imQTL_cpg_variant_df <- tmp_df %>% select(phenotype_id,variant_id) %>% mutate(combination=combination) %>% mutate(Dataset=Dataset)
    imQTL_cpg_variant_df <- rbind(imQTL_cpg_variant_df,tmp_imQTL_cpg_variant_df)
  }
}
# removing eosinophils
parsed_imQTL_cpg_variant_df <- imQTL_cpg_variant_df %>% separate(combination,remove=T,sep="_",into=c("F1","F2","Tissue","celltype_precursor","F3","F4","F5")) %>% separate(celltype_precursor, sep= "\\.", into=c("Celltype","F6")) %>% select(Dataset, Tissue, Celltype, phenotype_id, variant_id) %>% filter(!grepl("Eosino",Celltype))
# retaining only wb, colon, and lung
parsed_imQTL_cpg_variant_df <- parsed_imQTL_cpg_variant_df %>% filter(Tissue %in% c("wb","colon","lung"))
# creating vectors of imQTL elements
imqtl_cpg_GTEx_colon <- unique((parsed_imQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="colon"))$phenotype_id)
imqtl_variant_GTEx_colon <- unique((parsed_imQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="colon"))$variant_id)
imqtl_cpg_GTEx_lung <- unique((parsed_imQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="lung"))$phenotype_id)
imqtl_variant_GTEx_lung <- unique((parsed_imQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="lung"))$variant_id)
imqtl_cpg_HEALS_wb <- unique((parsed_imQTL_cpg_variant_df %>% filter(Dataset=="HEALS") %>% filter(Tissue=="wb"))$phenotype_id)
imqtl_variant_HEALS_wb <- unique((parsed_imQTL_cpg_variant_df %>% filter(Dataset=="HEALS") %>% filter(Tissue=="wb"))$variant_id)

#######################################
# extracting cpg and variant IDs for mQTLs
mQTL_cpg_variant_df <- data.frame()
for (Dataset in c("GTEx","HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/mQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl.txt.gz",list.files())]
  for (file_name in file_list) {
    print(file_name)
    tmp_df <- fread(file_name)
    tmp_df <- tmp_df %>% mutate(pval_adj_bh=p.adjust(pval_beta,method="fdr")) %>% filter(pval_adj_bh < 0.05)
    
    # assembling a DF containing all the cpg and variant IDs
    tmp_mQTL_cpg_variant_df <- tmp_df %>% select(phenotype_id,variant_id) %>% mutate(file_name=file_name) %>% mutate(Tissue=gsub(".cis_qtl.txt.gz","",file_name)) %>% mutate(Tissue=gsub("tensorQTL_mQTL_","",Tissue)) %>% mutate(Dataset=Dataset)
    mQTL_cpg_variant_df <- rbind(mQTL_cpg_variant_df,tmp_mQTL_cpg_variant_df)
  }
}
# creating vectors of mQTL elements
mqtl_cpg_GTEx_colon <- unique((mQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="colon"))$phenotype_id)
mqtl_variant_GTEx_colon <- unique((mQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="colon"))$variant_id)
mqtl_cpg_GTEx_lung <- unique((mQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="lung"))$phenotype_id)
mqtl_variant_GTEx_lung <- unique((mQTL_cpg_variant_df %>% filter(Dataset=="GTEx") %>% filter(Tissue=="lung"))$variant_id)
mqtl_cpg_HEALS_wb <- unique((mQTL_cpg_variant_df %>% filter(Dataset=="HEALS") %>% filter(Tissue=="wb"))$phenotype_id)
mqtl_variant_HEALS_wb <- unique((mQTL_cpg_variant_df %>% filter(Dataset=="HEALS") %>% filter(Tissue=="wb"))$variant_id)





# identify all bed types to examine enrichment for TBFS
background_bed_list <- list.files("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed")
pooled_estimates_list <- vector(mode="list",length=length(background_bed_list))
for (bed_index in 1:length(background_bed_list)) {
  background_bed=background_bed_list[bed_index]
  core_name <- gsub("^all_|\\.bed$", "", background_bed)
  print(core_name)
  imqtl_elements <- get(paste0("imqtl_", core_name))
  mqtl_elements <- get(paste0("mqtl_", core_name))
  background_elements <- unique((fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed/",background_bed)))$V4)
  
  # displaying total number of each element type
  print(paste("Processing:", core_name))
  num_imqtl_elements<-length(unique(imqtl_elements))
  num_mqtl_elements<-length(unique(mqtl_elements))
  num_background_elements<-length(unique(background_elements))
  
  # importing background tfbs
  import_background_tfbs_df <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/bedtools_output/all_",core_name,".bed-encRegTfbsClustered.bed"))
  # identifying tfbs for each qtl type
  imqtl_tfbs_df <- import_background_tfbs_df %>% filter(V4%in%imqtl_elements)
  mqtl_tfbs_df <- import_background_tfbs_df %>% filter(V4%in%mqtl_elements)
  background_imqtl_tfbs_df <- import_background_tfbs_df %>% filter(V8%in%imqtl_tfbs_df$V8) 
  background_mqtl_tfbs_df <- import_background_tfbs_df %>% filter(V8%in%mqtl_tfbs_df$V8) 
  
  #################################
  # Enrichment of imqtls for tfbs #
  #################################
  # obtaining counts of tfbs for those identified in imqtls
  imqtl_tfbs_count <- data.frame(table(imqtl_tfbs_df$V8))%>%rename(count_imqtl_tfbs=Freq)
  # obtaining counts of background tfbs for those identified in imqtls
  background_imqtl_tfbs_count <- data.frame(table(background_imqtl_tfbs_df$V8))%>%rename(count_background_tfbs=Freq)
  # joining tables and adding counts 
  join_imqtl_tfbs_count <- inner_join(imqtl_tfbs_count,background_imqtl_tfbs_count) %>% mutate(count_imqtl_elements=num_imqtl_elements,count_background_elements=num_background_elements)
  # computing enrichment scores
  join_imqtl_tfbs_count <- join_imqtl_tfbs_count %>%
    rowwise() %>%
    mutate(
      a = ifelse(count_imqtl_tfbs == 0, 0.5, as.numeric(count_imqtl_tfbs)),
      b = ifelse((count_imqtl_elements - count_imqtl_tfbs) == 0, 0.5, as.numeric(count_imqtl_elements - count_imqtl_tfbs)),
      c = ifelse(count_background_tfbs == 0, 0.5, as.numeric(count_background_tfbs)),
      d = ifelse((count_background_elements - count_background_tfbs) == 0, 0.5, as.numeric(count_background_elements - count_background_tfbs)),
      beta_imqtl = log((a * d) / (b * c)),
      beta_imqtl = ifelse(is.infinite(beta_imqtl), NA, beta_imqtl),
      se_imqtl = sqrt(1 / a + 1 / b + 1 / c + 1 / d),
      p_imqtl = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
    ) %>%
    ungroup() %>% data.frame()
  
  #################################
  # Enrichment of mqtls for tfbs #
  #################################
  # obtaining counts of tfbs for those identified in mqtls
  mqtl_tfbs_count <- data.frame(table(mqtl_tfbs_df$V8))%>%rename(count_mqtl_tfbs=Freq)
  # obtaining counts of background tfbs for those identified in mqtls
  background_mqtl_tfbs_count <- data.frame(table(background_mqtl_tfbs_df$V8))%>%rename(count_background_tfbs=Freq)
  # joining tables and adding counts 
  join_mqtl_tfbs_count <- inner_join(mqtl_tfbs_count,background_mqtl_tfbs_count) %>% mutate(count_mqtl_elements=num_mqtl_elements,count_background_elements=num_background_elements)
  # computing enrichment scores
  join_mqtl_tfbs_count <- join_mqtl_tfbs_count %>%
    rowwise() %>%
    mutate(
      a = ifelse(count_mqtl_tfbs == 0, 0.5, as.numeric(count_mqtl_tfbs)),
      b = ifelse((count_mqtl_elements - count_mqtl_tfbs) == 0, 0.5, as.numeric(count_mqtl_elements - count_mqtl_tfbs)),
      c = ifelse(count_background_tfbs == 0, 0.5, as.numeric(count_background_tfbs)),
      d = ifelse((count_background_elements - count_background_tfbs) == 0, 0.5, as.numeric(count_background_elements - count_background_tfbs)),
      beta_mqtl = log((a * d) / (b * c)),
      beta_mqtl = ifelse(is.infinite(beta_mqtl), NA, beta_mqtl),
      se_mqtl = sqrt(1 / a + 1 / b + 1 / c + 1 / d),
      p_mqtl = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
    ) %>%
    ungroup() %>% data.frame()
  
  #########################################
  # Pooling estimates from imqtl and mqtl #
  #########################################
  pooled_estimates <- inner_join(join_imqtl_tfbs_count%>%select(Var1,beta_imqtl,se_imqtl,p_imqtl),join_mqtl_tfbs_count%>%select(Var1,beta_mqtl,se_mqtl,p_mqtl),by=c("Var1"))
  
  names(pooled_estimates_list)[bed_index]<-core_name
  pooled_estimates_list[[bed_index]] <- pooled_estimates
}



############################
# PERFORMING META-ANALYSIS #
############################
library(dplyr)
library(tidyr)
library(purrr)
library(meta)

# ---- Step 1: Define Suffixes to Standardize ----
allowed_suffixes <- c("beta_imqtl", "se_imqtl", "p_imqtl", "beta_mqtl", "se_mqtl", "p_mqtl")

# ---- Step 2: Cleanly Rename Columns with Tissue Prefix ----
rename_with_tissue_prefix <- function(df, tissue_prefix) {
  new_names <- sapply(colnames(df), function(col) {
    if (col == "Var1") {
      return("Var1")
    }
    for (sfx in allowed_suffixes) {
      if (endsWith(col, sfx)) {
        return(paste0(tissue_prefix, "_", sfx))
      }
    }
    return(paste0(tissue_prefix, "_", col))  # fallback
  })
  colnames(df) <- new_names
  return(df)
}

# ---- Step 3: Meta-analysis Function ----
run_meta <- function(beta_vec, se_vec) {
  keep <- which(!is.na(beta_vec) & !is.na(se_vec))
  n_keep <- length(keep)
  
  if (n_keep >= 2) {
    res <- metagen(
      TE = beta_vec[keep],
      seTE = se_vec[keep],
      sm = "SMD",
      common = FALSE,
      random = TRUE
    )
    return(c(beta = res$TE.random, se = res$seTE.random, p = res$pval.random))
  } else if (n_keep == 1) {
    b <- beta_vec[keep]
    s <- se_vec[keep]
    p <- 2 * pnorm(abs(b / s), lower.tail = FALSE)
    return(c(beta = b, se = s, p = p))
  } else {
    return(c(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
}

# ---- Step 4: Define Tissue Prefix Mapping ----
tissue_map <- list(
  cpg_GTEx_colon     = "colon",
  cpg_GTEx_lung      = "lung",
  cpg_HEALS_wb       = "wb",
  variant_GTEx_colon = "colon",
  variant_GTEx_lung  = "lung",
  variant_HEALS_wb   = "wb"
)

# ---- Step 5: Create New List With Renamed Columns ----
renamed_estimates_list <- imap(pooled_estimates_list, function(df, name) {
  prefix <- tissue_map[[name]]
  rename_with_tissue_prefix(df, prefix)
})

# ---- Step 6: Join CpG and Variant Data Separately ----
cpg_joined <- reduce(
  renamed_estimates_list[startsWith(names(renamed_estimates_list), "cpg_")],
  full_join, by = "Var1"
)

variant_joined <- reduce(
  renamed_estimates_list[startsWith(names(renamed_estimates_list), "variant_")],
  full_join, by = "Var1"
)

# ---- Step 7: Row-wise Meta-analysis Application ----
perform_meta <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      meta_imqtl = list(run_meta(
        beta_vec = unlist(c_across(matches("_beta_imqtl$"))),
        se_vec   = unlist(c_across(matches("_se_imqtl$")))
      )),
      meta_mqtl = list(run_meta(
        beta_vec = unlist(c_across(matches("_beta_mqtl$"))),
        se_vec   = unlist(c_across(matches("_se_mqtl$")))
      ))
    ) %>%
    ungroup() %>%
    mutate(
      meta_beta_imqtl = sapply(meta_imqtl, function(x) x[1]),
      meta_se_imqtl   = sapply(meta_imqtl, function(x) x[2]),
      meta_p_imqtl    = sapply(meta_imqtl, function(x) x[3]),
      meta_beta_mqtl  = sapply(meta_mqtl, function(x) x[1]),
      meta_se_mqtl    = sapply(meta_mqtl, function(x) x[2]),
      meta_p_mqtl     = sapply(meta_mqtl, function(x) x[3])
    ) %>%
    select(-meta_imqtl, -meta_mqtl)
}

# ---- Step 8: Run Meta-analysis ----
meta_cpg_df     <- perform_meta(cpg_joined) %>% as.data.frame()
meta_variant_df <- perform_meta(variant_joined) %>% as.data.frame()

# ---- Step 9: Extract Clean Final Results ----
meta_cpg_df_final <- meta_cpg_df %>%
  select(Var1, starts_with("meta_")) %>%
  as.data.frame()
meta_variant_df_final <- meta_variant_df %>%
  select(Var1, starts_with("meta_")) %>%
  as.data.frame()

# collating all results prior to plotting
all_tfbs_results <- rbind(
  meta_cpg_df %>% mutate(element="CpG sites"),
  meta_variant_df %>% mutate(element="Variants")
) %>% as.data.frame() %>% rename(annotation=Var1)

# reshape into long format
long_tfbs_results <- all_tfbs_results %>%
  pivot_longer(
    cols = matches("^(colon|lung|wb|meta)_(beta|se|p)_(imqtl|mqtl)$"),
    names_to = c("tissue", ".value"),
    names_pattern = "(colon|lung|wb|meta)_(beta_imqtl|se_imqtl|p_imqtl|beta_mqtl|se_mqtl|p_mqtl)"
  ) %>%
  relocate(tissue, .after = annotation) %>% mutate(group="TFBS") %>% as.data.frame()

# making sure only to track enrichment -- thus, filtering for positive effect sizes in imQTLs
long_tfbs_results <- long_tfbs_results %>% filter(beta_imqtl > 0)
head(long_tfbs_results)
write.table(long_tfbs_results,file="/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/TFBS.tsv",quote=F,sep="\t",row.names=F,col.names=T)
# writing out enriched TFBS
write.table(long_tfbs_results%>%filter(element=="Variants")%>%filter(p_imqtl<0.05),file="/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/SIG_TFBS_VARIANTS.tsv",quote=F,sep="\t",row.names=F,col.names=T)


##################################################################
##### TFBS ENRICHMENT FOREST PLOTS FOR BOTH CPG AND VARIANTS #####
##################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(data.table)
library(patchwork)

# Output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/plots/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to generate a centered placeholder plot
make_placeholder_plot <- function(el_type) {
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    annotate("text", x = 0, y = 0,
             label = paste0("No statistically significant\nenriched TFBS among ", el_type),
             size = 3.5, hjust = 0.5, vjust = 0.5) +
    theme_void()
}

for (tissue_name in c("colon", "lung", "wb", "meta")) {
  
  plot_df <- long_tfbs_results %>%
    filter(tissue == tissue_name,
           !is.na(beta_imqtl),
           !is.na(se_imqtl),
           !is.na(beta_mqtl),
           !is.na(se_mqtl)) %>%
    mutate(
      lower_imqtl = beta_imqtl - 1.96 * se_imqtl,
      upper_imqtl = beta_imqtl + 1.96 * se_imqtl,
      lower_mqtl = beta_mqtl - 1.96 * se_mqtl,
      upper_mqtl = beta_mqtl + 1.96 * se_mqtl
    ) %>%
    pivot_longer(
      cols = c(beta_imqtl, beta_mqtl, lower_imqtl, lower_mqtl, upper_imqtl, upper_mqtl),
      names_to = c(".value", "type"),
      names_pattern = "(beta|lower|upper)_(imqtl|mqtl)"
    ) %>%
    mutate(
      type = recode(type, imqtl = "imQTL", mqtl = "mQTL"),
      annotation = as.factor(annotation)
    )
  
  element_plots <- lapply(c("CpG sites", "Variants"), function(el_type) {
    
    filtered_df <- long_tfbs_results %>%
      filter(tissue == tissue_name,
             element == el_type,
             !is.na(p_imqtl))
    
    bonf_threshold <- 0.05 / nrow(filtered_df)
    
    top_annotations <- filtered_df %>%
      filter(p_imqtl < bonf_threshold) %>%
      arrange(p_imqtl) %>%
      distinct(annotation, .keep_all = TRUE) %>%
      slice_head(n = 15) %>%
      select(annotation, group, p_imqtl_ranked = p_imqtl)
    
    if (nrow(top_annotations) == 0) {
      return(make_placeholder_plot(el_type))
    }
    
    df_el <- plot_df %>%
      filter(element == el_type, annotation %in% top_annotations$annotation) %>%
      left_join(top_annotations, by = c("annotation", "group")) %>%
      group_by(annotation) %>%
      mutate(beta_imqtl_value = beta[type == "imQTL"][1]) %>%
      ungroup() %>%
      mutate(annotation = fct_reorder(annotation, beta_imqtl_value, .desc = TRUE)) %>%
      droplevels() %>%
      group_by(group) %>%
      mutate(annotation = fct_rev(annotation)) %>%
      ungroup()
    
    if (nrow(df_el) == 0) {
      return(make_placeholder_plot(el_type))
    }
    
    # Stripe shading
    strip_rects <- df_el %>%
      distinct(group, annotation) %>%
      group_by(group) %>%
      mutate(
        y = rev(row_number()),
        ymin = y - 0.5,
        ymax = y + 0.5,
        shade = rep(c("white", "gray95"), length.out = n())
      ) %>%
      ungroup()
    
    ggplot(df_el, aes(x = beta, y = annotation, color = type)) +
      geom_rect(data = strip_rects, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = shade),
                inherit.aes = FALSE, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(position = position_dodge(width = 0.7)) +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0,
                     position = position_dodge(width = 0.7)) +
      facet_grid(rows = vars(group), scales = "free_y", space = "free_y", switch = "y", drop = FALSE) +
      scale_color_manual(values = c("imQTL" = "red", "mQTL" = "blue")) +
      scale_fill_identity() +
      labs(
        title = el_type,
        x = "Log Odds Ratio",
        y = NULL,
        color = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "bold"),
        axis.text.y = element_text(size = 11),
        panel.spacing = unit(1, "lines")
      )
  })
  
  combined_plot <- element_plots[[1]] + element_plots[[2]] + plot_layout(ncol = 2)
  
  ggsave(
    filename = paste0(output_dir, paste0("tfbs_",tolower(gsub(" ", "_", tissue_name)), ".png")),
    plot = combined_plot,
    width = 6, height = 5, dpi = 300
  )
}

##########################################
##### FOREST PLOTS FOR ONLY VARIANTS #####
##########################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(data.table)
library(patchwork)

# Output directory
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/plots/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Placeholder plot
make_placeholder_plot <- function(el_type) {
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    annotate("text", x = 0, y = 0,
             label = paste0("No statistically significant\nenriched TFBS among ", el_type),
             size = 3.5, hjust = 0.5, vjust = 0.5) +
    theme_classic()
}

for (tissue_name in c("colon", "lung", "wb", "meta")) {
  
  plot_df <- long_tfbs_results %>%
    filter(tissue == tissue_name,
           !is.na(beta_imqtl),
           !is.na(se_imqtl),
           !is.na(beta_mqtl),
           !is.na(se_mqtl)) %>%
    mutate(
      lower_imqtl = beta_imqtl - 1.96 * se_imqtl,
      upper_imqtl = beta_imqtl + 1.96 * se_imqtl,
      lower_mqtl = beta_mqtl - 1.96 * se_mqtl,
      upper_mqtl = beta_mqtl + 1.96 * se_mqtl
    ) %>%
    pivot_longer(
      cols = c(beta_imqtl, beta_mqtl, lower_imqtl, lower_mqtl, upper_imqtl, upper_mqtl),
      names_to = c(".value", "type"),
      names_pattern = "(beta|lower|upper)_(imqtl|mqtl)"
    ) %>%
    mutate(
      type = recode(type, imqtl = "imQTL", mqtl = "mQTL"),
      annotation = as.factor(annotation)
    )
  
  # Only plot for "Variants"
  el_type <- "Variants"
  
  filtered_df <- long_tfbs_results %>%
    filter(tissue == tissue_name,
           element == el_type,
           !is.na(p_imqtl))
  
  bonf_threshold <- 0.05 / nrow(filtered_df)
  
  top_annotations <- filtered_df %>%
    filter(p_imqtl < bonf_threshold) %>%
    arrange(p_imqtl) %>%
    distinct(annotation, .keep_all = TRUE) %>%
    slice_head(n = 20) %>%
    select(annotation, group, p_imqtl_ranked = p_imqtl)
  
  if (nrow(top_annotations) == 0) {
    plot_obj <- make_placeholder_plot(el_type)
  } else {
    df_el <- plot_df %>%
      filter(element == el_type, annotation %in% top_annotations$annotation) %>%
      left_join(top_annotations, by = c("annotation", "group")) %>%
      group_by(annotation) %>%
      mutate(beta_imqtl_value = beta[type == "imQTL"][1]) %>%
      ungroup() %>%
      mutate(annotation = fct_reorder(annotation, beta_imqtl_value, .desc = TRUE)) %>%
      droplevels() %>%
      mutate(annotation = fct_rev(annotation))
    
    strip_rects <- df_el %>%
      distinct(annotation) %>%
      mutate(
        y = rev(row_number()),
        ymin = y - 0.5,
        ymax = y + 0.5,
        shade = rep(c("white", "gray95"), length.out = n())
      )
    
    plot_obj <- ggplot(df_el, aes(x = beta, y = annotation, color = type)) +
      geom_rect(data = strip_rects, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = shade),
                inherit.aes = FALSE, alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(position = position_dodge(width = 0.7)) +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0,
                     position = position_dodge(width = 0.7)) +
      scale_color_manual(values = c("imQTL" = "red", "mQTL" = "blue")) +
      scale_fill_identity() +
      labs(
        title = "",
        x = "Log Odds Ratio",
        y = NULL,
        color = NULL
      ) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.text.y = element_text(size = 13),
        panel.spacing = unit(1, "lines")
      )
  }
  
  ggsave(
    filename = paste0(output_dir, paste0("tfbs_variants_",tolower(gsub(" ", "_", tissue_name)), ".png")),
    plot = plot_obj,
    width = 4.75, height = 4.75, dpi = 300
  )
}
