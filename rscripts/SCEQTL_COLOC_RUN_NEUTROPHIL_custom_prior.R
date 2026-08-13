# Code to download sceQTL files for Neutrophils
# cd /gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/SCEQTL/blueprint
# wget https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv
# grep QTS000002 tabix_ftp_paths.tsv | grep "neutrophil"
# wget https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000026/QTD000026.all.tsv.gz



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
    tmp_imqtl_cpg_variant_df <- tmp_df %>% select(phenotype_id,variant_id) %>% mutate(combination=combination) %>% mutate(Dataset=Dataset)
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


# importing eqtls
eqtl <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/SCEQTL/blueprint/QTD000026.all.tsv.gz")
library(data.table)
setDT(eqtl)
tmp_eqtl <- eqtl[
  ,
  .(
    phenotype_id = molecular_trait_id,
    variant_id = chartr("_", ":", variant),
    maf = maf,
    sample_size = an / 2,
    pval_nominal = pvalue,
    slope = beta,
    slope_se = se
  )
]

# filtering imqtl table to just have those which are also eQTL variants
imqtl_filtered <- parsed_imqtl_cpg_variant_df %>% filter(variant_id %in% tmp_eqtl$variant_id) %>% mutate(file_variant_id = gsub("\\:","_",variant_id)) %>% mutate(extracted_coloc_input_path = paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_mqtl_data/",parsed_combination,"_",phenotype_id,"_",file_variant_id,".tsv"))




######################
## Performing coloc ##
######################

library(data.table)
library(coloc)

##################################################
# EXPECTED OBJECTS ALREADY LOADED:
# tmp_eqtl
# imqtl_filtered
##################################################

WINDOW <- 250000
MIN_SHARED <- 50

setDT(tmp_eqtl)
setDT(imqtl_filtered)

##################################################
# PREPROCESS EQTL TABLE ONCE
##################################################

tmp_eqtl[, c("eqtl_chr", "eqtl_pos", "a1", "a2") :=
           tstrsplit(variant_id, ":", fixed = TRUE)]

tmp_eqtl[, eqtl_pos := as.integer(eqtl_pos)]

tmp_eqtl[, eqtl_MAF := fifelse(maf > 0.5, 1 - maf, maf)]
tmp_eqtl[, eqtl_varbeta := slope_se^2]

setindex(tmp_eqtl, eqtl_chr, eqtl_pos)

imqtl_filtered[, pos := as.integer(pos)]

##################################################
# RUN COLOC
##################################################

results <- vector("list", nrow(imqtl_filtered))

for (i in seq_len(nrow(imqtl_filtered))) {
  
  message("Performing coloc on imQTL ", i, " out of ", nrow(imqtl_filtered))
  
  lead_cpg   <- imqtl_filtered$phenotype_id[i]
  lead_var   <- imqtl_filtered$variant_id[i]
  lead_chr   <- imqtl_filtered$chr[i]
  lead_pos   <- imqtl_filtered$pos[i]
  lead_combo <- imqtl_filtered$parsed_combination[i]
  imqtl_path <- imqtl_filtered$extracted_coloc_input_path[i]
  
  if (!file.exists(imqtl_path)) {
    warning("Missing imQTL file: ", imqtl_path)
    next
  }
  
  ##################################################
  # READ IMQTL REGION
  ##################################################
  
  tmp_imqtl <- fread(imqtl_path)
  
  tmp_imqtl <- tmp_imqtl[
    ,
    .(
      variant_id,
      imqtl_beta = slope,
      imqtl_varbeta = slope_se^2,
      imqtl_MAF = fifelse(af > 0.5, 1 - af, af),
      imqtl_N = ma_samples
    )
  ]
  
  tmp_imqtl <- na.omit(tmp_imqtl)
  
  if (!(lead_var %in% tmp_imqtl$variant_id)) {
    next
  }
  
  ##################################################
  # FILTER EQTLS NEAR LEAD IMQTL
  ##################################################
  
  lo <- lead_pos - WINDOW
  hi <- lead_pos + WINDOW
  
  eqtl_region <- tmp_eqtl[
    eqtl_chr == lead_chr &
      eqtl_pos >= lo &
      eqtl_pos <= hi
  ]
  
  if (nrow(eqtl_region) <= MIN_SHARED) {
    next
  }
  
  ##################################################
  # KEEP ONLY SHARED VARIANTS
  ##################################################
  
  eqtl_region <- eqtl_region[
    variant_id %in% tmp_imqtl$variant_id
  ]
  
  if (nrow(eqtl_region) <= MIN_SHARED) {
    next
  }
  
  ##################################################
  # RUN COLOC PER EQTL GENE
  ##################################################
  
  coloc_results <- list()
  genes <- unique(eqtl_region$phenotype_id)
  
  for (gene in genes) {
    
    eqtl_gene <- eqtl_region[phenotype_id == gene]
    
    shared_variants <- intersect(eqtl_gene$variant_id, tmp_imqtl$variant_id)
    
    in_eqtl  <- lead_var %in% eqtl_gene$variant_id
    in_imqtl <- lead_var %in% tmp_imqtl$variant_id
    
    if (length(shared_variants) <= MIN_SHARED || !in_eqtl || !in_imqtl) {
      next
    }
    
    dataset1 <- eqtl_gene[
      variant_id %in% shared_variants,
      .(
        snp = variant_id,
        beta = slope,
        varbeta = eqtl_varbeta,
        MAF = eqtl_MAF,
        N = sample_size,
        type = "quant"
      )
    ]
    # removing any potential duplicate snps in the eqtl coloc input dataset
    dataset1 <- dataset1 %>%
      distinct(snp, .keep_all = TRUE)
    
    dataset2 <- tmp_imqtl[
      variant_id %in% shared_variants,
      .(
        snp = variant_id,
        beta = imqtl_beta,
        varbeta = imqtl_varbeta,
        MAF = imqtl_MAF,
        N = imqtl_N,
        type = "quant"
      )
    ]
    
    dataset1 <- na.omit(dataset1)
    dataset2 <- na.omit(dataset2)
    
    shared_final <- intersect(dataset1$snp, dataset2$snp)
    
    dataset1 <- dataset1[snp %in% shared_final]
    dataset2 <- dataset2[snp %in% shared_final]
    
    if (length(shared_final) <= MIN_SHARED) {
      next
    }
    
    coloc_input1 <- list(
      snp = dataset1$snp,
      beta = dataset1$beta,
      varbeta = dataset1$varbeta,
      MAF = dataset1$MAF,
      N = median(dataset1$N, na.rm = TRUE),
      type = "quant"
    )
    
    coloc_input2 <- list(
      snp = dataset2$snp,
      beta = dataset2$beta,
      varbeta = dataset2$varbeta,
      MAF = dataset2$MAF,
      N = median(dataset2$N, na.rm = TRUE),
      type = "quant"
    )
    
    coloc_out <- coloc.abf(coloc_input1, coloc_input2, p1 = 0.0005, p2 = 0.0045, p12 = 0.0005)
    
    coloc_results[[gene]] <- data.table(
      parsed_combination = lead_combo,
      cpg_id = lead_cpg,
      imqtl_variant_id = lead_var,
      chr = lead_chr,
      pos = lead_pos,
      eqtl_phenotype_id = gene,
      nsnps = as.numeric(coloc_out$summary["nsnps"]),
      PP.H0.abf = as.numeric(coloc_out$summary["PP.H0.abf"]),
      PP.H1.abf = as.numeric(coloc_out$summary["PP.H1.abf"]),
      PP.H2.abf = as.numeric(coloc_out$summary["PP.H2.abf"]),
      PP.H3.abf = as.numeric(coloc_out$summary["PP.H3.abf"]),
      PP.H4.abf = as.numeric(coloc_out$summary["PP.H4.abf"])
    )
  }
  
  coloc_summary <- rbindlist(coloc_results, fill = TRUE)
  
  if (nrow(coloc_summary) == 0) {
    next
  }
  
  results[[i]] <- coloc_summary
}

##################################################
# SAVE RESULTS
##################################################

results_df <- rbindlist(results, fill = TRUE)

setorder(results_df, -PP.H4.abf)

# adding sceqtl cell type column
results_df <- results_df %>% mutate(sceqtl_ct="Neutro")

fwrite(
  results_df,
  "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/custom_prior/coloc_output_Neutro.tsv",
  sep = "\t"
)

