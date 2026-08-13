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

################################################
# importing imqtl, mqtl, and background elements and determine how many intersect each annotation #
################################################
background_bed_list <- list.files("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed")
enrichment_df <- data.frame()
for (background_bed in background_bed_list) {
  core_name <- gsub("^all_|\\.bed$", "", background_bed)
  imqtl_elements <- get(paste0("imqtl_", core_name))
  mqtl_elements <- get(paste0("mqtl_", core_name))
  background_elements <- unique((fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed/",background_bed)))$V4)
  
  # displaying total number of each element type
  print(paste("Processing:", core_name))
  num_imqtl_elements<-length(unique(imqtl_elements))
  num_mqtl_elements<-length(unique(mqtl_elements))
  num_background_elements<-length(unique(background_elements))
  
  # iterating through target elements and identifying overlap
  target_file_list <- list.files("/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/bedtools_output/")[grepl(background_bed,list.files("/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/bedtools_output/"))]
  for (current_target in target_file_list) {
    #current_target=target_file_list[1]
    current_target_elements<-(fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/bedtools_output/",current_target)))$V4
    
    # identifying number of intersecting elements with the current annotation
    num_intersect_imqtl_elements<-length(unique(intersect(imqtl_elements,current_target_elements)))
    num_intersect_mqtl_elements<-length(unique(intersect(mqtl_elements,current_target_elements)))
    num_intersect_background_elements<-length(unique(intersect(background_elements,current_target_elements)))
    
    # extracting the annotation string
    annotation <- str_split(current_target, "-", simplify = TRUE)[2]
    
    # storing results in df
    tmp_enrichment_df <- data.frame(core_name,annotation,num_imqtl_elements,num_intersect_imqtl_elements,num_mqtl_elements,num_intersect_mqtl_elements,num_background_elements,num_intersect_background_elements)
    enrichment_df <- rbind(enrichment_df, tmp_enrichment_df)
  }
}

############################################## 
########## PARSING ENRICHMENT TABLE ########## 
############################################## 
# removing mis-matched tissues in the chromatin segmentation state analysis
parsed_enrichment_df <- enrichment_df %>%
  filter(!(
    (grepl("colon", core_name) & grepl("wb", annotation)) |
      (grepl("wb", core_name) & grepl("colon", annotation)) |
      (grepl("colon", core_name) & grepl("lung", annotation)) |
      (grepl("lung", core_name) & grepl("colon", annotation)) |
      (grepl("wb", core_name) & grepl("lung", annotation)) |
      (grepl("lung", core_name) & grepl("wb", annotation))
  ))

# further parsing enrichment table
parsed_enrichment_df <- parsed_enrichment_df %>% separate(core_name,into=c("element","dataset","tissue"),sep="_") %>% mutate(annotation=gsub("colon_","",gsub("wb_","",gsub("lung_","",gsub(".bed","",annotation))))) %>%
  mutate(annotation = case_when(
    annotation == "encode5_distalenhancers" ~ "Distal Enhancer",
    annotation == "encode5_promoterlikeregions" ~ "Promoter",
    annotation == "encode5_insulators" ~ "Insulator",
    annotation == "encode5_proximalenhancers" ~ "Proximal Enhancer",
    annotation == "wgEncodeGencodeBasicV26_TSS1500" ~ "TSS1500",
    annotation == "wgEncodeGencodeBasicV26_TSS200" ~ "TSS200",
    annotation == "wgEncodeGencodeBasicV26_cds" ~ "CDS",
    annotation == "wgEncodeGencodeBasicV26_intron" ~ "Intron",
    annotation == "wgEncodeGencodeBasicV26_utr3" ~ "3' UTR",
    annotation == "wgEncodeGencodeBasicV26_utr5" ~ "5' UTR",
    annotation == "wgEncodeRegDnaseClustered" ~ "DNaseI HS",
    annotation == "encRegTfbsClustered" ~ "TFBS",
    annotation == "cpgIsland" ~ "CGI",
    annotation == "cpgShore" ~ "CGI/Shores",
    annotation == "cpgShelf" ~ "CGI/Shores/Shelves",
    TRUE ~ annotation  # keep other annotations unchanged
  ))

# making groupings and parsing chromatin states
chromatin_states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF_Rpts", "9_Het",
                      "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
chromatin_clean <- str_remove(chromatin_states, "^[0-9]+_")
cpg_regions <- c("CGI", "CGI/Shores", "CGI/Shores/Shelves")
reg_regions <- c("TFBS", "Distal Enhancer", "Insulator", "Promoter", "Proximal Enhancer", "DNaseI HS")
gene_contexts <- c("TSS1500", "TSS200", "CDS", "Intron", "3' UTR", "5' UTR")

parsed_enrichment_df <- parsed_enrichment_df %>%
  mutate(annotation = if_else(annotation %in% chromatin_states, str_remove(annotation, "^[0-9]+_"), annotation)) %>% 
  mutate(group = case_when(
    annotation %in% chromatin_clean ~ "Chromatin States",
    annotation %in% cpg_regions ~ "CpG Island Regions",
    annotation %in% reg_regions ~ "Regulatory Regions",
    annotation %in% gene_contexts ~ "Gene Contexts",
    TRUE ~ "Other")) %>%
  mutate(annotation=ifelse(annotation=="ZNF_Rpts","ZNF/Rpts",annotation)) %>%
  mutate(element = case_when(
    element == "cpg" ~ "CpG sites",
    element == "variant" ~ "Variants",
    TRUE ~ element  # keep other annotations unchanged
  )) %>%
  mutate(tissue = case_when(
    tissue == "colon" ~ "Colon",
    tissue == "lung" ~ "Lung",
    tissue == "wb" ~ "Whole Blood",
    TRUE ~ tissue  # keep other annotations unchanged
  ))








##########################################################
# COMPUTING ENRICHMENT OF EACH ANNOTATION IN EACH TISSUE #
##########################################################
############## UTILIZING FISHER EXACT TESTS ##############
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

# storing as a df indicating it has enrichment reults by tissue
bytissue_enrichment_df <- parsed_enrichment_df
write.table(x = bytissue_enrichment_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/bytissue_enrichment_df.tsv",row.names=F,quote=F,col.names=T,sep="\t")


####################################################
# META-ANALYZING ENRICHMENT RESULTS ACROSS TISSUES #
####################################################
library(dplyr)
library(tidyr)
library(meta)

# Step 1: Clean tissue names to avoid space-related column naming issues
parsed_enrichment_df_removespacetissue <- parsed_enrichment_df %>%
  mutate(tissue = make.names(tissue))  # e.g., "Whole Blood" -> "Whole.Blood"

# Step 2: Pivot to wide format
wide_df <- parsed_enrichment_df_removespacetissue %>%
  filter(tissue %in% c("Colon", "Lung", "Whole.Blood")) %>%
  select(element, annotation, group, tissue,
         beta_imqtl, se_imqtl, beta_mqtl, se_mqtl) %>%
  pivot_wider(
    names_from = tissue,
    values_from = c(beta_imqtl, se_imqtl, beta_mqtl, se_mqtl)
  )

# Step 3: Updated meta-analysis function with single-tissue fallback
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
    # Use single tissue values
    b <- beta_vec[keep]
    s <- se_vec[keep]
    p <- 2 * pnorm(abs(b / s), lower.tail = FALSE)
    return(c(beta = b, se = s, p = p))
    
  } else {
    return(c(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
}

# Step 4: Apply meta-analysis row-wise
meta_results <- wide_df %>%
  rowwise() %>%
  mutate(
    meta_beta_imqtl = run_meta(c_across(starts_with("beta_imqtl_")),
                               c_across(starts_with("se_imqtl_")))[1],
    meta_se_imqtl   = run_meta(c_across(starts_with("beta_imqtl_")),
                               c_across(starts_with("se_imqtl_")))[2],
    meta_p_imqtl    = run_meta(c_across(starts_with("beta_imqtl_")),
                               c_across(starts_with("se_imqtl_")))[3],
    
    meta_beta_mqtl = run_meta(c_across(starts_with("beta_mqtl_")),
                              c_across(starts_with("se_mqtl_")))[1],
    meta_se_mqtl   = run_meta(c_across(starts_with("beta_mqtl_")),
                              c_across(starts_with("se_mqtl_")))[2],
    meta_p_mqtl    = run_meta(c_across(starts_with("beta_mqtl_")),
                              c_across(starts_with("se_mqtl_")))[3]
  ) %>%
  ungroup()

# Step 5: Return clean data.frame with just the meta-analysis results
meta_enrichment_df <- meta_results %>%
  select(element, annotation, group,
         meta_beta_imqtl, meta_se_imqtl, meta_p_imqtl,
         meta_beta_mqtl, meta_se_mqtl, meta_p_mqtl) %>%
  as.data.frame()

write.table(x = meta_enrichment_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/ENRICHMENT_ANALYSIS/tables/meta_enrichment_df.tsv",row.names=F,quote=F,col.names=T,sep="\t")

