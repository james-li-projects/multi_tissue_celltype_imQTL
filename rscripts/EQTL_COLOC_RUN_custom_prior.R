library(dplyr)
library(data.table)
library(coloc)
library(tidyr)
library(purrr)
library(tibble)

############################################
# IMPORTING ARGUMENTS FOR IMPORTING IMQTLS #
############################################
# importing list of imqtls to process
args <- commandArgs(trailingOnly = TRUE)
imported_imqtl_query_name <- args[1]
# imported_imqtl_query_name <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_list/query_1.tsv"
imported_imqtl_query <- fread(imported_imqtl_query_name)
imqtl_file_list <- (imported_imqtl_query %>% mutate(filelist=paste0(
  parsed_combination,"_",
  phenotype_id,"_",
  gsub("\\:","_",variant_id),".tsv"
)) %>% select(filelist) )$filelist
total_num_imqtl <- length(imqtl_file_list)

###############################################
# INITIALIZING TISSUE LIST AND OUTPUT OBJECTS #
###############################################
# initializing tissue list
tissue_list <- c("colon","lung","ovary","prostate","wb")
tissue_column_names <- as.vector(rbind(
  paste0(tissue_list, "_ensg"),
  paste0(tissue_list, ".PP.H0.abf"),
  paste0(tissue_list, ".PP.H1.abf"),
  paste0(tissue_list, ".PP.H2.abf"),
  paste0(tissue_list, ".PP.H3.abf"),
  paste0(tissue_list, ".PP.H4.abf")
))

# initializing matrix to store coloc results
result_matrix <- matrix(NA, nrow = length(imqtl_file_list), ncol = 6 * length(tissue_list) + 5)
colnames(result_matrix) <- c("combination","cpg_id","variant_id","chr","pos", tissue_column_names)

########################
# INITALIZE LOOP IMQTL #
########################
# iterating through all imQTLs 
for (i in 1:length(imqtl_file_list)) {
  
  # printing out evaluated imqtl
  print(paste("Performing coloc on imqtl",i,"out of",total_num_imqtl))
  
  # importing file argument for current imQTL file
  tmp_imqtl_file <- imqtl_file_list[i]
  
  # parsing out relevant strings from imQTL filename
  tmp_cleaned <- sub("\\.tsv$", "", tmp_imqtl_file)
  parts <- strsplit(tmp_cleaned, "_")[[1]]
  tmp_imqtl_combination <- paste(parts[1:2], collapse = "_")
  tmp_imqtl_chr_int <- as.integer(gsub("chr","",parts[4]))
  tmp_imqtl_pos_int <- as.integer(parts[5])
  tmp_imqtl_cpg_id <- parts[3]
  tmp_imqtl_variant_id <- paste(parts[4:7], collapse = ":")
  
  #######################################
  # IMPORTING AND PROCESSING EACH IMQTL #
  #######################################
  # importing imqtl DF
  tmp_imqtl <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_data/",tmp_imqtl_file))
  tmp_imqtl <- tmp_imqtl %>% mutate(imqtl_N=round(ma_count/(2*af)),imqtl_MAF=ifelse(af>0.5,1-af,af),imqtl_varbeta=b_gi_se^2,imqtl_beta=b_gi,imqtl_type="quant") %>% select(variant_id,imqtl_beta,imqtl_varbeta,imqtl_MAF,imqtl_N,phenotype_id,imqtl_type)
  tmp_imqtl <- na.omit(tmp_imqtl)
  # storing imqtl annotations
  result_matrix[i,1] <- tmp_imqtl_combination
  result_matrix[i,2] <- tmp_imqtl_cpg_id
  result_matrix[i,3] <- tmp_imqtl_variant_id
  result_matrix[i,4] <- tmp_imqtl_chr_int
  result_matrix[i,5] <- tmp_imqtl_pos_int
  
  ###################################################
  # IMPORTING AND PROCESSING EQTLS FROM EACH TISSUE #
  ###################################################
  # initializing tissue index
  for (m in 1:length(tissue_list)) {
    tmp_tissue=tissue_list[m]
    tmp_eqtl <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GTEx_eQTL_v8/output/tsv/tensorQTL_eQTL_",tmp_tissue,".cis_qtl_pairs.chr",tmp_imqtl_chr_int,".tsv"))
    tmp_eqtl[, c("chr", "pos", "a2", "a1") := tstrsplit(variant_id, ":", fixed=TRUE)]
    tmp_eqtl[, pos := as.integer(pos)]
    tmp_eqtl <- tmp_eqtl[abs(pos - tmp_imqtl_pos_int) <= 250000]
    
    # if no intersecting eqtls, then skip coloc analysis for this tissue and go to the next eqtl tissue
    if (nrow(tmp_eqtl>50)) {
      tmp_eqtl <- tmp_eqtl %>%
        mutate(
          sample_size = round(ma_count / (2 * af))
        )
      tmp_eqtl <- na.omit(tmp_eqtl)
      
      # dividing up eqtls by gene (phenotype_id)
      split_phenotype_ids <- unique(tmp_eqtl$phenotype_id)
      eqtl_list_named <- list()
      for (pid in split_phenotype_ids) {
        eqtl_list_named[[as.character(pid)]] <- tmp_eqtl %>% filter(phenotype_id == pid)
      }
      
      # Initialize empty list to store results
      coloc_results <- list()
      
      # Get all phenotype names (names of the list)
      phenotype_ids <- names(eqtl_list_named)
      
      # Begin loop
      for (k in seq_along(phenotype_ids)) {
        phenotype <- phenotype_ids[k]
        tmp_eqtl_subset <- eqtl_list_named[[k]]
        
        # Identify shared variants between this phenotype and imQTL
        shared_variants <- intersect(tmp_eqtl_subset$variant_id, tmp_imqtl$variant_id)
        
        if (length(shared_variants) > 50) {
          # Dataset 1: eQTL from tmp_eqtl
          dataset1 <- tmp_eqtl_subset %>%
            filter(variant_id %in% shared_variants) %>%
            mutate(
              snp = variant_id,
              beta = slope,
              varbeta = slope_se^2,
              MAF = ifelse(af>0.5,1-af,af),
              N = sample_size,
              type = "quant"
            ) %>%
            select(snp, beta, varbeta, MAF, N, type)
          
          # Dataset 2: imQTL
          dataset2 <- tmp_imqtl %>%
            filter(variant_id %in% shared_variants) %>%
            mutate(
              snp = variant_id,
              beta = imqtl_beta,
              varbeta = imqtl_varbeta,
              MAF = imqtl_MAF,
              N = imqtl_N,
              type = "quant"
            ) %>%
            select(snp, beta, varbeta, MAF, N, type)
          
          # Build coloc input lists
          coloc_input1 <- list(
            snp     = dataset1$snp,
            beta    = dataset1$beta,
            varbeta = dataset1$varbeta,
            MAF     = dataset1$MAF,
            N       = median(dataset1$N),
            type    = unique(dataset1$type)
          )
          
          coloc_input2 <- list(
            snp     = dataset2$snp,
            beta    = dataset2$beta,
            varbeta = dataset2$varbeta,
            MAF     = dataset2$MAF,
            N       = median(dataset2$N),
            type    = unique(dataset2$type)
          )
          
          # Run coloc
          coloc_out <- coloc.abf(coloc_input1, coloc_input2, p1 = 0.0005, p2 = 0.0045, p12 = 0.0005)
          PP.H0.abf <- coloc_out$summary["PP.H0.abf"]
          PP.H1.abf <- coloc_out$summary["PP.H1.abf"]
          PP.H2.abf <- coloc_out$summary["PP.H2.abf"]
          PP.H3.abf <- coloc_out$summary["PP.H3.abf"]
          PP.H4.abf <- coloc_out$summary["PP.H4.abf"]
          
          # Save result
          coloc_results[[phenotype]] <- tibble(
            phenotype_id = phenotype, 
            PP.H0.abf = PP.H0.abf,
            PP.H1.abf = PP.H1.abf,
            PP.H2.abf = PP.H2.abf,
            PP.H3.abf = PP.H3.abf,
            PP.H4.abf = PP.H4.abf)
        }
      }
      
      # Combine all non-null results into a single dataframe
      coloc_summary <- bind_rows(coloc_results)
      
      if (nrow(coloc_summary) > 0) {
        # obtaining gene name and all PP.abf for the gene with the largest PP.H4.abf
        highest_coloc <- coloc_summary %>% arrange(desc(PP.H4.abf)) %>% head(1)
        highest_ensg <- as.character(highest_coloc$phenotype_id)
        highest_PP.H0.abf <- as.numeric(highest_coloc$PP.H0.abf)
        highest_PP.H1.abf <- as.numeric(highest_coloc$PP.H1.abf)
        highest_PP.H2.abf <- as.numeric(highest_coloc$PP.H2.abf)
        highest_PP.H3.abf <- as.numeric(highest_coloc$PP.H3.abf)
        highest_PP.H4.abf <- as.numeric(highest_coloc$PP.H4.abf)
        
        # storing results
        result_matrix[i,(m-1)*6+1+5] <- highest_ensg
        result_matrix[i,(m-1)*6+2+5] <- highest_PP.H0.abf
        result_matrix[i,(m-1)*6+3+5] <- highest_PP.H1.abf
        result_matrix[i,(m-1)*6+4+5] <- highest_PP.H2.abf
        result_matrix[i,(m-1)*6+5+5] <- highest_PP.H3.abf
        result_matrix[i,(m-1)*6+6+5] <- highest_PP.H4.abf
      }
    }
  }
}

# converting results matrix to a DF
results_df <- data.frame(data.frame(result_matrix))

############################
# OUTPUTTING COLOC RESULTS #
############################
# generating output path
output_path <- gsub("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_list/","/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/coloc_output_custom_prior/",imported_imqtl_query_name)

# writing out coloc results as a tsv
write.table(results_df,file=output_path,quote=F,row.names=F,col.names=T,sep="\t")
