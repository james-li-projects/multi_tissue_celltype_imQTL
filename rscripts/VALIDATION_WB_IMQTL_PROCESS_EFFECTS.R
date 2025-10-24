library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(ggnewscale) 
library(dplyr)
library(tidyr)
library(purrr)

# set working directory
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/regional_plots")

# processing imQTLs across all tissue and cell type combinations
all_imqtl_results <- data.frame()
gtex_tissue_list <- c("wb")
heals_tissue_list <- c("wb")
for (dataset in c("GTEx","HEALS")) {
  if (dataset=="GTEx") {
    tissue_list <- gtex_tissue_list
    pvar<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix.pvar")
    afreq<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix_FREQ.afreq") %>% filter(ALT_FREQS>=0.1)
  } else if (dataset=="HEALS") {
    tissue_list <- heals_tissue_list
    pvar<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix.pvar")
    afreq<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix_FREQ.afreq") %>% filter(ALT_FREQS>=0.1)
  }
  
  # iterating across all tissues
  for (current_tissue in tissue_list) {
    # importing tissue DNA methylation bed file
    input_dnam_bed <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",dataset,"/",current_tissue,".bed"),sep="\t") 
    # importing tissue covariates files
    cov_df <- data.frame(t(read.table(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",dataset,"/processed_covariates_",current_tissue,".txt"),sep="\t"))) 
    if (dataset=="GTEx") {
      cov_df <- cov_df %>% mutate(Sample_Name=gsub("\\.","-",rownames(cov_df)))
    } else if (dataset=="HEALS") {
      cov_df <- cov_df %>% mutate(Sample_Name=gsub("^X","",rownames(cov_df)))
    }
    # write out tissue list of samples
    write.table(data.frame(0,cov_df$Sample_Name),file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list",quote=F,row.names=F,col.names=F,sep="\t")
    
    # iterating across all cell type combinations
    setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",dataset,"/imQTL/top_assoc"))
    tmp_combination_list <- list.files()[grepl(current_tissue,list.files()) & grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
    tmp_combination_list <- gsub("tensorQTL_imQTL_","",gsub("\\.cis_qtl_top_assoc.txt.gz","",tmp_combination_list))
    for (tmp_combination in tmp_combination_list) {
      # importing CT interaction file
      ct_df <- read.table(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",dataset,"/processed_interactions_",tmp_combination,".txt"),sep="\t") %>% rename(Sample_Name=V1,CT=V2)
      
      # importing imQTL results [only HEALS imQTLs will be evaluated for validation
      tmp_imqtl_df <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc/tensorQTL_imQTL_",tmp_combination,".cis_qtl_top_assoc.txt.gz")) %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni")) %>% filter(pval_adj_bonf<0.05)
      if (nrow(tmp_imqtl_df) < 1) {
        print("No imQTLs detected at Bonferroni<0.05")
      } else {
        print(paste("Current tissue/CT combination:",tmp_combination))
        
        # iterating across all imQTLs (Bonferroni<0.05) detected for the current tissue cell type combination
        for (b in 1:nrow(tmp_imqtl_df)) {
          
          # initializing target cpg and lead variant
          target_cpg = tmp_imqtl_df$phenotype_id[b]
          lead_variant_raw = tmp_imqtl_df$variant_id[b]

          # making sure the cpg and snp are both present, otherwise skipping this record
          tmp_size_cpg <- nrow(input_dnam_bed %>% filter(phenotype_id==target_cpg) %>% head(1))
          tmp_size_variant <- nrow(afreq %>% filter(ID==lead_variant_raw))
          indicator_progress=tmp_size_cpg+tmp_size_variant
          if (indicator_progress==2) {
          
          ###############################
          # reformatting lead variant with periods for later matching
          lead_variant_period<-gsub("\\:",".",lead_variant_raw)
          
          # importing all variants that should be utilized in regional plot
          num_chr<-as.numeric(gsub("chr","",(stringr::str_split(lead_variant_raw,pattern="\\:"))[[1]][1]))
          num_pos<-as.numeric(gsub("chr","",(stringr::str_split(lead_variant_raw,pattern="\\:"))[[1]][2]))
          window_variants <- pvar %>% filter(`#CHROM`==paste0("chr",num_chr)) %>% filter(POS>(num_pos-10000),POS<(num_pos+10000))
          write.table(window_variants %>% select(ID),file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list",quote=F,row.names=F,col.names=F,sep="\t")
          
          # filering DNAm bed file 
          dnam_bed <- input_dnam_bed %>% filter(phenotype_id==target_cpg) %>% head(1)
          t_dnam_bed <- data.frame(t(dnam_bed))
          t_dnam_bed <- t_dnam_bed %>% slice(-c(1:4)) 
          cpg_df <- t_dnam_bed %>% rename(mval=t.dnam_bed.) %>% mutate(Sample_Name=rownames(t_dnam_bed)) %>% select(Sample_Name,everything())
          
          # importing variant data
          system(paste0("module load plink/2.0; plink2 -pfile /gpfs/data/pierce-lab/james.li/imQTL/data/",dataset,"/genetic_data/processed_genetic_data_chrprefix --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list --keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list --maf 0.05 --export Av --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window"))
          traw <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window.traw") %>% separate(SNP,into=c("T1","T2","T3","T4"),remove=F) %>% mutate(across(11:ncol(.), ~ ifelse(COUNTED != T4, 2 - ., .))) %>% select(-c(paste0(rep("T","4"),1:4)))

          # obtaining df genetic variants
          variant_df <- 
            data.frame(t(traw))
          newcolnames <- variant_df["SNP",]
          variant_df <- variant_df %>%
            rename_with(~ as.character( newcolnames ), everything())
          # Filter out the rows to be removed
          rows_to_remove <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
          variant_df <- variant_df %>%
            filter(!rownames(.) %in% rows_to_remove)
          # adding sample name column
          variant_df <- variant_df %>% mutate(Sample_Name = gsub("^0_","",rownames(.))) %>% select(Sample_Name,everything())
          
          ########################
          # Join the data frames by "Sample_Name" using inner_join
          reg_df <- cpg_df %>%
            inner_join(ct_df, by = "Sample_Name") %>%
            inner_join(cov_df, by = "Sample_Name") %>%
            inner_join(variant_df, by = "Sample_Name") %>% select(-Sample_Name)
          colnames(reg_df) <- gsub("\\:",".",colnames(reg_df))
          
          # Convert all predictors to numeric
          convert_to_numeric <- function(df) {
            df[] <- lapply(df, function(col) {
              if (is.numeric(col)) {
                return(col)
              } else {
                # Convert non-numeric columns to numeric
                return(as.numeric(col))
              }
            })
            return(df)
          }
          # Ensure all predictors are numeric
          reg_df <- convert_to_numeric(reg_df)
          
          # Identify column indices of reg_df that start with "chr"
          chr_indices <- grep("^chr", colnames(reg_df))
          
          # Identify the index of the 'lead_variant'
          lead_variant_index <- which(gsub("\\:",".",colnames(reg_df)) == lead_variant_period)
          
          # Get indices of columns starting with "chr" but not equal to 'lead_variant_raw'
          evaluated_variant_indices <- setdiff(chr_indices, lead_variant_index)
          
          
          #######################################################
          # Function to run regression and extract coefficients #
          #######################################################
          run_regression <- function(df, outcome, predictors) {
            # Construct the formula
            formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
            
            # Generate model matrix to check rank deficiency
            model_matrix <- model.matrix(formula, data = df)
            model_rank <- qr(model_matrix)$rank  # Get rank of matrix
            expected_rank <- length(predictors) + 1  # Expected rank (including intercept)
            
            # If rank is lower than expected, return NA for all coefficients
            if (model_rank < expected_rank) {
              coef_names <- c("(Intercept)", predictors)  # Ensure correct variable names
              return(matrix(NA, nrow = length(coef_names), ncol = 4, 
                            dimnames = list(coef_names, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))))
            }
            
            # Run linear regression
            model <- lm(formula, data = df)
            
            # Extract coefficients from summary if no collinearity issue
            return(summary(model)$coefficients)
          }
          
          ################################################
          # INITIALIZING DATA FRAME TO STORE ALL RESULTS #
          ################################################
          results <- data.frame()
          
          ###################################################
          # OBTAINING MARGINAL, MAIN, AND INTERACTION TERMS #
          ###################################################
          # initializing list of baseline predictors and the outcome variable
          tmp_outcome <- "mval"
          tmp_baseline_predictors <- setdiff(colnames(reg_df[-c(chr_indices)]),"mval")
          for (current_variant_col_index in lead_variant_index) {
            # obtaining name of variant being analyzed
            current_variant_col_name <- colnames(reg_df)[current_variant_col_index]
            
            ############
            # Marginal #
            ############
            # obtaining list of predictors
            # tmp_predictors <- setdiff(c(current_variant_col_name,tmp_baseline_predictors),"CT")
            # obtaining regression value for current variant 
            #regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
            #tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Marginal Effect")
            #results<-rbind(results,tmp_results)
            
            ################################
            # Interaction and Main Effects #
            ################################
            # obtaining list of predictors
            tmp_interaction_term <- paste0(current_variant_col_name,"*CT")
            tmp_predictors <- c(current_variant_col_name,tmp_baseline_predictors,tmp_interaction_term)
            # obtaining regression value for current variant 
            regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
            
            ###############
            # Interaction #
            ###############
            tmp_results <- cbind(regression_results[gsub("\\*",":",tmp_interaction_term),],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Interaction Effect")
            results<-rbind(results,tmp_results)
            
            ########
            # Main #
            ########
            #tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Main Effect")
            #results<-rbind(results,tmp_results)
          }

          
          # last bug check, sometimes results have nothing
          if (nrow(results) > 0) {
            #################
            # FINAL PARSING #
            #################
            # assembling marginal/main/interaction plot DF
            results<-data.frame(results)
            rownames(results) <- NULL
            colnames(results) <- c("BETA","SE","T","P","variant_id","phenotype_id","class")
            results <- results %>% mutate(BETA=as.numeric(BETA),P=as.numeric(P),SE=as.numeric(SE),T=as.numeric(T)) 
            # adding a column for dataset
            results <- results %>% mutate(current_dataset=dataset) %>% rename(dataset=current_dataset)
            
            # storing results specifically for the interaction terms of imqtls
            tmp_imqtl_results <- results %>% filter(variant_id==lead_variant_period) %>% mutate(combination=tmp_combination) 
            all_imqtl_results <- rbind(all_imqtl_results,tmp_imqtl_results)
          }
          } # this block will progress if both cpg and variant are present in the dataset
          
        }
      }
    }
  }
}

# filtering out eosinophil imQTL
all_imqtl_results <- all_imqtl_results %>% filter(combination!="wb_Eosino")

# saving all parsed imqtl effect sizes
save(all_imqtl_results,file="/gpfs/data/pierce-lab/james.li/imQTL/output/validate_imqtl_wb_gtex/all_imqtl_results.RData")
