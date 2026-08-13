# Code to download sceQTL files for PBMCs
# wget--content-disposition https://zenodo.org/records/15097677/files/B_IN_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/B_Mem_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/CD4_ET_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/CD4_NC_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/CD4_SOX4_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/CD8_ET_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/CD8_NC_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/CD8_S100B_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/DC1_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/Mono_C_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/Mono_NC_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/NK_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/NK_R_JOBS_summary_statistics.txt.gz?download=1
# wget--content-disposition https://zenodo.org/records/15097677/files/Plasma_JOBS_summary_statistics.txt.gz?download=1

#######################################
# importing libraries
#######################################

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(coloc)

# importing european allele frequencies
eur.afreq<-fread("/gpfs/data/pierce-lab/james.li/EUR_1kg_hg38/eur.afreq")
setDT(eur.afreq)
eur.afreq <- eur.afreq[
  ,
  .(
    snp_hg38 = paste0(gsub(":", "_", ID, fixed = TRUE), "_b38"),
    MAF = pmin(ALT_FREQS, 1 - ALT_FREQS)
  )
]
setDT(eur.afreq)
eur.afreq <- unique(eur.afreq, by = "snp_hg38")

#######################################
# set working directory
#######################################

setwd("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/SCEQTL/40154479_Wang")

sceqtl_ct_list <- gsub("_JOBS_summary_statistics.txt.gz", "", list.files())

for (current_sceqtl_ct in sceqtl_ct_list) {
  
  #######################################
  # importing eqtl summary statistics for the current cell type
  #######################################
  
  eqtl <- fread(paste0(
    "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/SCEQTL/40154479_Wang/",
    current_sceqtl_ct,
    "_JOBS_summary_statistics.txt.gz"
  ))
  
  # joining eqtls to european allele frequencies from 1kg
  setDT(eqtl)
  setkey(eqtl, snp_hg38)
  setkey(eur.afreq, snp_hg38)
  eqtl <- eur.afreq[eqtl, nomatch = 0]
  
  # obtaining hg38 positions for each eQTL
  eqtl[, c("chrom", "pos_hg38", "a1", "a2") :=
         tstrsplit(snp_hg38, "_", fixed = TRUE, keep = 1:4)]
  
  eqtl[, `:=`(
    chrom = as.integer(chrom),
    pos_hg38 = as.integer(pos_hg38)
  )]
  
  #######################################
  # importing whole blood imQTLs
  #######################################
  
  Dataset <- c()
  combination <- c()
  imqtl_cpg_variant_df <- data.frame()
  
  for (Dataset in c("HEALS")) {
    
    setwd(paste0(
      "/gpfs/data/pierce-lab/james.li/imQTL/output/",
      Dataset,
      "/imQTL/top_assoc"
    ))
    
    file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz", list.files())]
    
    for (combination in file_list) {
      
      print(combination)
      
      tmp_df <- fread(combination)
      
      # filtering for Bonferroni-adjusted p-value < 0.05
      tmp_df <- tmp_df %>%
        mutate(pval_adj_bonf = p.adjust(pval_emt, method = "bonferroni")) %>%
        filter(pval_adj_bonf < 0.05)
      
      # assembling a DF containing all the CpG-variant pairs
      tmp_imqtl_cpg_variant_df <- tmp_df %>%
        select(phenotype_id, variant_id) %>%
        mutate(combination = combination) %>%
        mutate(Dataset = Dataset)
      
      imqtl_cpg_variant_df <- rbind(
        imqtl_cpg_variant_df,
        tmp_imqtl_cpg_variant_df
      )
    }
  }
  
  #######################################
  # parsing imQTL variant info
  #######################################
  
  parsed_imqtl_cpg_variant_df <- imqtl_cpg_variant_df %>%
    separate(
      variant_id,
      remove = FALSE,
      into = c("chr", "pos", "a2", "a1"),
      sep = "\\:"
    ) %>%
    mutate(chunk_num = as.integer(gsub("chr", "", chr))) %>%
    mutate(
      parsed_combination = gsub(
        "tensorQTL_imQTL_",
        "",
        gsub("\\.cis_qtl_top_assoc.txt.gz", "", combination)
      )
    ) %>%
    mutate(
      input_parquet_filename = paste0(
        "/gpfs/data/pierce-lab/james.li/imQTL/output/",
        Dataset,
        "/imQTL/full_assoc/tensorQTL_imQTL_",
        parsed_combination,
        ".chunk",
        chunk_num,
        ".cis_qtl_pairs.chr",
        chunk_num,
        ".parquet"
      )
    ) %>%
    select(-a2, -a1, -combination, -Dataset)
  
  # filtering out eosinophils
  parsed_imqtl_cpg_variant_df <- parsed_imqtl_cpg_variant_df %>%
    filter(parsed_combination != "wb_Eosino")
  
  #######################################
  # filtering eqtl table to just contain imQTL lead variant positions
  #######################################
  
  setDT(eqtl)
  setDT(parsed_imqtl_cpg_variant_df)
  
  # small lookup table
  lookup <- unique(parsed_imqtl_cpg_variant_df[, .(
    chr = chunk_num,
    pos_hg38 = pos
  )])
  
  # make sure types match
  eqtl[, chr := as.character(chr)]
  lookup[, chr := as.character(chr)]
  eqtl[, pos_hg38 := as.integer(pos_hg38)]
  lookup[, pos_hg38 := as.integer(pos_hg38)]
  
  # fast filter: keep eqtl rows whose chr + pos appear in lookup
  eqtl_filtered <- eqtl[lookup, on = .(chr, pos_hg38), nomatch = 0]
  
  # generating variant IDs to compare with imqtl results
  eqtl_filtered <- eqtl_filtered %>%
    mutate(
      ID1 = paste(paste0("chr", chrom), pos_hg38, a1, a2, sep = ":"),
      ID2 = paste(paste0("chr", chrom), pos_hg38, a2, a1, sep = ":")
    )
  
  #######################################
  # filtering imqtl table to just have those which are also eQTL variants
  #######################################
  
  imqtl_filtered <- parsed_imqtl_cpg_variant_df %>%
    filter(variant_id %in% eqtl_filtered$ID1 | variant_id %in% eqtl_filtered$ID2) %>%
    mutate(file_variant_id = gsub("\\:", "_", variant_id)) %>%
    mutate(
      extracted_coloc_input_path = paste0(
        "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_mqtl_data/",
        parsed_combination,
        "_",
        phenotype_id,
        "_",
        file_variant_id,
        ".tsv"
      )
    )
  
  #######################################
  # making eqtl data.table in line with coloc code
  #######################################
  
  setDT(eqtl)
  
  tmp_eqtl <- eqtl[
    MAF > 0 & MAF < 1,
    .(
      phenotype_id = gene,
      variant_id = paste0(
        "chr",
        gsub("_", ":", sub("_b38$", "", snp_hg38))
      ),
      pval_nominal = p_val,
      slope = beta.est,
      slope_se = beta.se,
      eqtl_MAF = fifelse(MAF > 0.5, 1 - MAF, MAF)
    )
  ]
  
  ####################
  # performing coloc #
  ####################
  
  EQTL_N <- 982
  WINDOW <- 250000
  MIN_SHARED <- 50
  
  setDT(tmp_eqtl)
  setDT(imqtl_filtered)
  
  tmp_eqtl[, c("eqtl_chr", "eqtl_pos", "a1", "a2") :=
             tstrsplit(variant_id, ":", fixed = TRUE)]
  
  tmp_eqtl[, eqtl_pos := as.integer(eqtl_pos)]
  
  setindex(tmp_eqtl, eqtl_chr, eqtl_pos)
  
  imqtl_filtered[, pos := as.integer(pos)]
  
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
    
    # Keep only variants that are also present in the imQTL region.
    # Note: this merge adds imqtl_MAF only for filtering/shared variant matching;
    # dataset1 below uses eqtl_MAF, not imqtl_MAF.
    eqtl_region <- merge(
      eqtl_region,
      tmp_imqtl[, .(variant_id, imqtl_MAF)],
      by = "variant_id",
      all = FALSE
    )
    
    if (nrow(eqtl_region) <= MIN_SHARED) {
      next
    }
    
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
          varbeta = slope_se^2,
          MAF = eqtl_MAF,
          N = EQTL_N,
          type = "quant"
        )
      ]
      
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
        N = EQTL_N,
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
  
  results_df <- rbindlist(results, fill = TRUE)
  
  setorder(results_df, -PP.H4.abf)
  
  # adding sceqtl cell type column
  results_df <- results_df %>%
    mutate(sceqtl_ct = current_sceqtl_ct)
  
  fwrite(
    results_df,
    paste0(
      "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/custom_prior/coloc_output_",
      current_sceqtl_ct,
      ".tsv"
    ),
    sep = "\t"
  )
}