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
gtex_tissue_list <- c("prostate","ovary","lung","colon")
heals_tissue_list <- c("wb")
for (dataset in c("GTEx","HEALS")) {
  if (dataset=="GTEx") {
    tissue_list <- gtex_tissue_list
    pvar<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix.pvar")
    
  } else if (dataset=="HEALS") {
    tissue_list <- heals_tissue_list
    pvar<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix.pvar")
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
      
      # importing imQTL results
      tmp_imqtl_df <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",dataset,"/imQTL/top_assoc/tensorQTL_imQTL_",tmp_combination,".cis_qtl_top_assoc.txt.gz")) %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni")) %>% filter(pval_adj_bonf<0.05)
      if (nrow(tmp_imqtl_df) < 1) {
        print("No imQTLs detected at Bonferroni<0.05")
      } else {
        print(paste("Current tissue/CT combination:",tmp_combination))
        # Create an output directory for the tissue cell type combination if it does not exist
        dir_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/regional_plots/comprehensive_plots/", tmp_combination)
        if (!dir.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
        }
        
        # iterating across all imQTLs (Bonferroni<0.05) detected for the current tissue cell type combination
        for (b in 1:nrow(tmp_imqtl_df)) {
          
          # initializing target cpg and lead variant
          target_cpg = tmp_imqtl_df$phenotype_id[b]
          lead_variant_raw = tmp_imqtl_df$variant_id[b]
          
          ###############################
          # reformatting lead variant with periods for later matching
          lead_variant_period<-gsub("\\:",".",lead_variant_raw)
          
          # importing all variants that should be utilized in regional plot
          num_chr<-as.numeric(gsub("chr","",(stringr::str_split(lead_variant_raw,pattern="\\:"))[[1]][1]))
          num_pos<-as.numeric(gsub("chr","",(stringr::str_split(lead_variant_raw,pattern="\\:"))[[1]][2]))
          window_variants <- pvar %>% filter(`#CHROM`==paste0("chr",num_chr)) %>% filter(POS>(num_pos-250000),POS<(num_pos+250000))
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
          for (current_variant_col_index in chr_indices) {
            # obtaining name of variant being analyzed
            current_variant_col_name <- colnames(reg_df)[current_variant_col_index]
            
            ############
            # Marginal #
            ############
            # obtaining list of predictors
            tmp_predictors <- setdiff(c(current_variant_col_name,tmp_baseline_predictors),"CT")
            # obtaining regression value for current variant 
            regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
            tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Marginal Effect")
            results<-rbind(results,tmp_results)
            
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
            tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Main Effect")
            results<-rbind(results,tmp_results)
          }
          
          ###############################################################
          # OBTAINING CONDITIONAL MARGINAL, MAIN, AND INTERACTION TERMS #
          ###############################################################
          # obtaining marginal effects
          tmp_outcome <- "mval"
          tmp_baseline_predictors <- setdiff(colnames(reg_df[-c(chr_indices)]),"mval")
          for (current_variant_col_index in evaluated_variant_indices) {
            # obtaining name of variant being analyzed
            current_variant_col_name <- colnames(reg_df)[current_variant_col_index]
            
            ############
            # Marginal #
            ############
            # obtaining list of predictors including the lead variant
            tmp_predictors <- c(setdiff(c(current_variant_col_name,tmp_baseline_predictors),"CT"),lead_variant_period)
            # obtaining regression value for current variant 
            regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
            tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Marginal Effect (conditional)")
            results<-rbind(results,tmp_results)
            
            ################################
            # Interaction and Main Effects #
            ################################
            # obtaining list of predictors including interaction terms with both the lead and target variants
            tmp_interaction_term_current_variant <- paste0(current_variant_col_name,"*CT")
            tmp_interaction_term_lead_variant <- paste0(lead_variant_period,"*CT")
            tmp_predictors <- c(current_variant_col_name,lead_variant_period,tmp_baseline_predictors,tmp_interaction_term_current_variant,tmp_interaction_term_lead_variant)
            # obtaining regression value for current variant 
            regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
            
            ###############
            # Interaction #
            ###############
            tmp_results <- cbind(regression_results[gsub("\\*",":",tmp_interaction_term_current_variant),],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Interaction Effect (conditional)")
            results<-rbind(results,tmp_results)
            
            ########
            # Main #
            ########
            tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Main Effect (conditional)")
            results<-rbind(results,tmp_results)
          }
          
          
          ####################################
          # Add r2 correlations between SNPs #
          ####################################
          # computing LD
          max_snp_window<-nrow(window_variants)
          system(paste0("module load plink/1.9; plink -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/",dataset,"/genetic_data/processed_genetic_data_chrprefix_bfile --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list --keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list --maf 0.05 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_LD"))
          
          # import pairwise correlations
          pairwise_cor <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_LD.ld")
          
          # Filter the rows where either SNP_A or SNP_B equals lead_variant_raw
          filtered_df <- pairwise_cor[pairwise_cor$SNP_A == lead_variant_raw | pairwise_cor$SNP_B == lead_variant_raw, ]
          
          # Create a new column 'OTHER_VARIANT'
          filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == lead_variant_raw, filtered_df$SNP_B, filtered_df$SNP_A)
          
          # finalizing r2 df
          r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(variant_id=OTHER_VARIANT,r2=R2)
          r2_df <- rbind(r2_df,data.frame(variant_id=lead_variant_raw,r2=1))
          r2_df <- r2_df %>% mutate(variant_id=gsub("\\:",".",variant_id))
          
          ##########################
          # PREPARING FOR PLOTTING #
          ##########################
          # assembling marginal/main/interaction plot DF
          regional_plot_df<-data.frame(results)
          rownames(regional_plot_df) <- NULL
          colnames(regional_plot_df) <- c("BETA","SE","T","P","variant_id","phenotype_id","class")
          # making an indicator variable for whether an observation is about the lead variant
          regional_plot_df <- regional_plot_df %>% mutate(isleadvariant=ifelse(variant_id==lead_variant_period,1,0))
          # preparing and parsing the table
          regional_plot_df <- regional_plot_df %>% mutate(BETA=as.numeric(BETA),P=as.numeric(P),SE=as.numeric(SE),T=as.numeric(T),Sign = ifelse(BETA > 0, "Positive", ifelse(BETA < 0, "Negative", NA))) %>% tidyr::separate(variant_id, into = c("chr_char","pos_char","a2","a1"),sep="\\.",remove=F) %>% mutate(chr_num=as.numeric(gsub("chr","",chr_char)),pos_num=as.numeric(pos_char)) %>% mutate(
            class=factor(regional_plot_df$class, 
                         levels = c("Marginal Effect", 
                                    "Marginal Effect (conditional)", 
                                    "Main Effect", 
                                    "Main Effect (conditional)", 
                                    "Interaction Effect", 
                                    "Interaction Effect (conditional)"))
          )
          
          # assembling LD plot DF
          regional_plot_df_2 <- regional_plot_df %>% filter(class %in% c("Interaction Effect","Interaction Effect (conditional)")) %>% mutate(
            class=ifelse(class=="Interaction Effect","LD plot - Interaction Effect","LD plot - Interaction Effect (conditional)")
          ) %>% mutate(class=factor(class,levels=c("LD plot - Interaction Effect","LD plot - Interaction Effect (conditional)"))) %>% inner_join(r2_df,by=c("variant_id"))
          
          ###################################
          # PLOTTING THE COMPREHENSIVE PLOT #
          ###################################
          library(ggplot2)
          library(ggnewscale)
          library(scales)
          
          # Determine the global min and max of Position for consistent x-axis
          x_min <- min(window_variants$POS)
          x_max <- max(window_variants$POS)
          
          # Create the plot
          plot <- ggplot() +
            
            # First set of points: Sign (Top 6 Facets)
            geom_point(data = regional_plot_df, 
                       aes(x = pos_num / 1e6, 
                           y = -log10(P), 
                           color = Sign, 
                           shape = factor(isleadvariant),
                           size = ifelse(isleadvariant == 1, 4, 3)),  # Increase diamond size
            ) +
            scale_color_manual(name = "Sign", values = c("Positive" = "blue", "Negative" = "red"),na.translate = FALSE) + 
            scale_shape_manual(values = c("0" = 16, "1" = 18)) +  # Circle for 0, Diamond for 1
            scale_size_identity() +  # Ensures size mapping is used correctly
            
            # Remove legend for shape (Lead Variant)
            guides(shape = "none") +  
            
            # Add a new color scale BEFORE adding the next geom_point()
            new_scale_color() + 
            
            # Second set of points: r2 (Bottom 2 Facets)
            geom_point(data = regional_plot_df_2, 
                       aes(x = pos_num / 1e6, 
                           y = -log10(P), 
                           color = r2, 
                           shape = factor(isleadvariant),
                           size = ifelse(isleadvariant == 1, 4, 3)),  # Increase diamond size
            ) +
            scale_color_gradient(name = expression(r^2), low = "yellow", high = "purple") +
            scale_shape_manual(values = c("0" = 16, "1" = 18)) +  # Reapply shape mapping
            scale_size_identity() +  # Ensures size mapping is used correctly
            
            # Remove legend for shape (Lead Variant) again to ensure it's gone
            guides(shape = "none") +
            
            # Facet by class in 2 columns
            facet_wrap(~ class, ncol = 2) +
            
            # Axis and plot labels
            labs(y = "-log10(P)") +
            
            # X-axis settings
            scale_x_continuous(
              name = paste0("Position on Chromosome ", num_chr, " (Mb)"),
              labels = comma_format(),
              breaks = c(x_min, x_max) / 1e6
            ) +
            
            # Theme settings
            theme_classic() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              text = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              panel.spacing = unit(0.5, "cm"),
              axis.text.x = element_text(hjust = c(0, 1))
            )
          
          # Save the plot
          ggsave(
            filename = paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/regional_plots/comprehensive_plots/",tmp_combination,"/regional_plot_",tmp_combination,"_", 
                              gsub("\\:", ".", lead_variant_period), ".png"),
            plot = plot,
            width = 10,
            height = 8
          )

          # storing results specifically for the imqtl (marginal/main/interaction effect sizes)
          tmp_imqtl_results <- results %>% filter(variant_id==lead_variant_period) %>% mutate(combination=tmp_combination) 
          all_imqtl_results <- rbind(all_imqtl_results,tmp_imqtl_results)
          
          # saving the regional_plot_df
          save(regional_plot_df,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/all_regional_plot_df/regional_plot_df_",tmp_combination,"_",lead_variant_period,"_",target_cpg,".RData"))
        }
      }
    }
  }
}

# filtering out eosinophil imQTL
all_imqtl_results <- all_imqtl_results %>% filter(combination!="wb_Eosino")

# saving all parsed imqtl effect sizes
save(all_imqtl_results,file="/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/all_imqtl_results.RData")
# load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/all_imqtl_results.RData")

# parsing this output DF 
parsed_imqtl_results <- all_imqtl_results %>% separate(combination,into=c("tissue","celltype"),remove=F,sep="_") %>% rename(`Regression Term`=class)


########################################
########################################
########### PARSING RESULTS ############
########################################
########################################
# Parsing long format to rename columns 
long_parsed_imqtl_results <- parsed_imqtl_results %>% rename(SE=`Std..Error`,P=`Pr...t..`)

# Convert to wide format
wide_parsed_imqtl_results <- long_parsed_imqtl_results %>%
  pivot_wider(
    id_cols = c(variant_id, phenotype_id, tissue, celltype, combination),
    names_from = `Regression Term`,
    values_from = c(Estimate, SE, t.value, P),
    names_sep = "_"
  ) %>%
  rename_with(~ gsub(" ", "", .x, fixed = TRUE)) %>%  # Remove whitespace
  mutate(across(everything(), ~ ifelse(map_lgl(.x, is.list), map(.x, unlist), .x)))  # Unlist list-columns

# Ensure all columns are atomic vectors
wide_parsed_imqtl_results <- wide_parsed_imqtl_results %>%
  mutate(across(where(is.list), ~ sapply(.x, function(x) ifelse(length(x) == 1, x, NA))))

# filtering the transformed dataframe to remove NA interaction term values
wide_parsed_imqtl_results <- wide_parsed_imqtl_results %>% filter(!is.na(P_InteractionEffect))

# Categorize the rows
wide_parsed_imqtl_results <- wide_parsed_imqtl_results %>%
  mutate(Category = case_when(
    # using a stringent 1e-5 threshold to categorize imQTLs as having a detectable mQTL marginal effect
    P_MarginalEffect > 1e-5 ~ "Unknown",
    (Estimate_MarginalEffect > 0 & Estimate_InteractionEffect > 0) |
      (Estimate_MarginalEffect < 0 & Estimate_InteractionEffect < 0) ~ "Same",
    (Estimate_MarginalEffect > 0 & Estimate_InteractionEffect < 0) |
      (Estimate_MarginalEffect < 0 & Estimate_InteractionEffect > 0) ~ "Different"
  ))

# Modify the dataframe based on the given conditions
consistent_wide_parsed_imqtl_results <- data.frame(wide_parsed_imqtl_results) %>%
  mutate(
    # reversing sign of interaction and main effects if marginal effect is negative
    Estimate_InteractionEffect = ifelse(Estimate_MarginalEffect < 0, -Estimate_InteractionEffect, Estimate_InteractionEffect),
    Estimate_MainEffect = ifelse(Estimate_MarginalEffect < 0, -Estimate_MainEffect, Estimate_MainEffect),
    t.value_InteractionEffect = ifelse(Estimate_MarginalEffect < 0, -t.value_InteractionEffect, t.value_InteractionEffect),
    t.value_MainEffect = ifelse(Estimate_MarginalEffect < 0, -t.value_MainEffect, t.value_MainEffect),
    # making marginal effects positive
    Estimate_MarginalEffect = abs(Estimate_MarginalEffect),  # Make it positive
    t.value_MarginalEffect = abs(t.value_MarginalEffect)  # Make it positive
  ) %>% mutate(
    `Interaction_Term_Sign`=ifelse(Estimate_InteractionEffect > 0,"(+) Interaction Sign","(-) Interaction Sign")
  )

summary(consistent_wide_parsed_imqtl_results$Estimate_MainEffect)
summary(consistent_wide_parsed_imqtl_results$Estimate_MarginalEffect)

# Reshape the dataframe into long format
consistent_long_parsed_imqtl_results <- data.frame(consistent_wide_parsed_imqtl_results %>%
                                                     pivot_longer(
                                                       cols = starts_with("Estimate_") | starts_with("SE_") | 
                                                         starts_with("t.value_") | starts_with("P_"),
                                                       names_to = c(".value", "reg_term"),
                                                       names_sep = "_"
                                                     )) %>% rename(`Regression Term` = reg_term) %>% mutate(
                                                       `Regression Term`=gsub("Effect"," Effect",`Regression Term`)
                                                     )
# Annotating whether each CpG was mapped to an mQTL
# Define tissue list
tissue_list <- c("colon", "kidney", "lung", "ovary","wb")
# Initialize an empty list to store results
result_list <- list()
# Loop through each tissue and process the corresponding file
for (current_tissue in tissue_list) {
  if (current_tissue!="wb") {
    directory <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/mQTL/top_assoc"
  } else if (current_tissue=="wb") {
    directory <- "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/mQTL/top_assoc"
  } else {}
  # Construct file path
  file_path <- file.path(directory, paste0("tensorQTL_mQTL_", current_tissue, ".cis_qtl.txt.gz"))
  # Read the data using fread
  dt <- fread(file_path)
  # Filter rows where pval_beta < 0.05
  dt_filtered <- dt[pval_beta < 0.05]
  # Extract unique phenotype_id values
  unique_phenotypes <- unique(dt_filtered$phenotype_id)
  # Store results in a data.table
  result_list[[current_tissue]] <- data.table(tissue = current_tissue, phenotype_id = unique_phenotypes)
}
# Combine results into a single data.table
final_dt <- rbindlist(result_list)
# Print or save the final result
print(final_dt)

# Add ismqtl column based on matching tissue and phenotype_id
consistent_long_parsed_imqtl_results <- consistent_long_parsed_imqtl_results %>%
  mutate(ismqtl = ifelse(paste(tissue, phenotype_id) %in% 
                           paste(final_dt$tissue, final_dt$phenotype_id), 1, 0))

# Display the first few rows to verify
head(consistent_long_parsed_imqtl_results)

# Add ismqtl column based on matching tissue and phenotype_id
consistent_wide_parsed_imqtl_results <- consistent_wide_parsed_imqtl_results %>%
  mutate(ismqtl = ifelse(paste(tissue, phenotype_id) %in% 
                           paste(final_dt$tissue, final_dt$phenotype_id), 1, 0))

# Display the first few rows to verify
head(consistent_wide_parsed_imqtl_results)


#######################################
############ PLOT BOXPLOTS ############
#######################################
# Create the box plot for imQTLs that are 1) yes mQTL/not mQTL, 2) consistent or not consistent interaction & marginal effects
library(dplyr)
library(ggplot2)
library(tidyr)

consistency_list <- c("(-) Interaction Sign", "(+) Interaction Sign")

for (j in 1:length(consistency_list)) {
  tmp_consistency <- consistency_list[j]
  
  # Define celltype map
  celltype_map <- c(
    "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
    "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
    "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
    "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
    "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
    "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
    "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
  )
  
  # Define tissue map
  tissue_map <- c(
    "colon" = "Colon", "lung" = "Lung", "prostate" = "Prostate",
    "ovary" = "Ovary", "wb" = "Whole Blood"
  )
  
  # Clean and prepare the data
  boxplot_df <- consistent_long_parsed_imqtl_results %>%
    filter(Interaction_Term_Sign == tmp_consistency) %>%
    mutate(
      clean_celltype = recode(celltype, !!!celltype_map),
      clean_tissue = recode(tissue, !!!tissue_map),
      tissue_celltype = paste0(clean_tissue, " - ", clean_celltype)
    ) %>%
    mutate(tissue_celltype = factor(tissue_celltype, levels = unique(tissue_celltype)))
  
  # Compute medians and group-wise max
  med_df <- boxplot_df %>%
    group_by(tissue_celltype, `Regression Term`) %>%
    summarize(
      med = median(Estimate, na.rm = TRUE),
      ymax = max(Estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_label = ymax + 0.10 * (max(boxplot_df$Estimate, na.rm = TRUE) - min(boxplot_df$Estimate, na.rm = TRUE)))
  
  # Plot
  p <- ggplot(boxplot_df, aes(x = tissue_celltype, y = Estimate, col = `Regression Term`)) +
    geom_boxplot() +
    geom_text(
      data = med_df,
      aes(x = tissue_celltype, y = y_label, label = sprintf("%.2f", med), color = `Regression Term`),
      position = position_dodge(width = 0.75),
      size = 2.5,
      inherit.aes = FALSE
    ) +
    labs(y = "Effect size estimate", x = "") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      plot.margin = margin(c(10, 10, 10, 50))
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray")
  
  # Save plot with wide width
  output_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/BOXPLOT_COMBINATIONS_consistency", j - 1, ".png")
  ggsave(output_path, plot = p, width = 20, height = 6, dpi = 600)
  cat("Box plot saved to:", output_path, "\n")
  
  # ---- Identify where "Marginal Effect" > "Main Effect" ----
  median_matrix <- med_df %>%
    select(tissue_celltype, `Regression Term`, med) %>%
    pivot_wider(names_from = `Regression Term`, values_from = med)
  
  if (all(c("Marginal Effect", "Main Effect") %in% colnames(median_matrix))) {
    higher_marginal <- median_matrix %>%
      filter(`Marginal Effect` > `Main Effect`) %>%
      pull(tissue_celltype)
    
    if (length(higher_marginal) > 0) {
      cat("Tissue-celltype combinations where Marginal Effect > Main Effect for consistency index", j - 1, ":\n")
      print(higher_marginal)
    } else {
      cat("No combinations where Marginal Effect > Main Effect for consistency index", j - 1, ".\n")
    }
  } else {
    cat("Warning: Could not find both 'Marginal Effect' and 'Main Effect' for consistency index", j - 1, "\n")
  }
}


#######################################
############ PLOT BARPLOTS ############
#######################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Define celltype map
celltype_map <- c(
  "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
  "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
  "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
  "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
  "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
  "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
  "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
)

# Clean up tissue names
consistent_wide_parsed_imqtl_results_label <- consistent_wide_parsed_imqtl_results %>%
  mutate(tissue = recode(tissue,
                         colon = "Colon",
                         lung = "Lung",
                         prostate = "Prostate",
                         ovary = "Ovary",
                         wb = "Whole Blood"))

# Extract celltype from combination and map it
consistent_wide_parsed_imqtl_results_label <- consistent_wide_parsed_imqtl_results_label %>%
  mutate(celltype = str_extract(combination, "[^_]+$")) %>%  # Extract the last segment of combination
  mutate(celltype_display = recode(celltype, !!!celltype_map))  # Map using celltype_map

# Set tissue order
tissue_order <- rev(c("Prostate", "Ovary", "Lung", "Colon", "Whole Blood"))

# Count per tissue+combo+category
category_counts <- consistent_wide_parsed_imqtl_results_label %>%
  count(tissue, combination, celltype_display, Category) %>%
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) %>%
  rename(Same = `Same`, Different = `Different`, Unknown = `Unknown`) %>%
  mutate(Label = paste0("(", Same, ", ", Different, ", ", Unknown, ")"),
         Total = Same + Different + Unknown)

# Order by tissue then Total
category_counts <- category_counts %>%
  mutate(tissue = factor(tissue, levels = tissue_order)) %>%
  arrange(tissue, desc(Total)) %>%
  mutate(combination = factor(combination, levels = unique(combination)))  # preserve order

# Long format for stacked bars
category_counts_long <- category_counts %>%
  pivot_longer(cols = c(Same, Different, Unknown), names_to = "Category", values_to = "n")

category_counts_long$Category <- factor(category_counts_long$Category, 
                                        levels = c("Same", "Different", "Unknown"))

# Plot
max_x <- max(category_counts$Total, na.rm = TRUE)
max_label_length <- max(nchar(category_counts$Label), na.rm = TRUE)

p <- ggplot(category_counts_long, aes(x = n, y = combination, fill = Category)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(data = category_counts, 
            aes(x = Total + max_x * 0.1 + max_label_length * 0.3, y = combination, label = Label), 
            inherit.aes = FALSE, size = 4, hjust = 0, color = "black") + 
  scale_fill_manual(name = "Directional Consistency",
                    values = c("Same" = "#00AB66", "Different" = "#FF8C00", "Unknown" = "gray")) +
  scale_y_discrete(labels = category_counts$celltype_display) +  # use mapped labels
  theme_classic() +
  labs(x = "Number of mapped imQTLs", y = NULL, fill = "Category") +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  ) +
  expand_limits(x = max_x * 1.5 + max_label_length * 2)

# Save
ggsave("/gpfs/data/pierce-lab/james.li/imQTL/output/directional_consistency/STACKED_BARPLOT.png", 
       plot = p, width = 4.5, height = 7, dpi = 300)

# checking consistency between these parsed results and tensorQTL output
combination_list <- unique(parsed_imqtl_results$combination)
combination_list <- combination_list[grepl("wb",combination_list)]
for (current_combination in combination_list) {
  print(current_combination)
  tensor_output<-fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc/tensorQTL_imQTL_",current_combination,".cis_qtl_top_assoc.txt.gz")) %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni")) %>% select(phenotype_id,variant_id,pval_gi,b_gi)
  
  r_output <- parsed_imqtl_results %>% filter(combination==current_combination) %>% filter(`Regression Term`=="Interaction Effect") %>% rename(pval_gi=`Pr...t..`,b_gi=Estimate) %>% mutate(variant_id=gsub("\\.",":",variant_id)) %>% select(variant_id,phenotype_id,pval_gi,b_gi)
  
  compare_df <- inner_join(tensor_output,r_output,by=c("phenotype_id","variant_id"))
  
  fit<-summary(lm(compare_df$b_gi.x~compare_df$b_gi.y))
  print(fit$r.squared)
  print(head(compare_df)%>%select(-phenotype_id))
}

####################################
# writing out parsed results 
wide_parsed_imqtl <- data.frame(wide_parsed_imqtl_results) %>% mutate(variant_id=gsub("\\.",":",variant_id))
save(wide_parsed_imqtl,file="/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")
write.table(wide_parsed_imqtl,file="/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.txt",quote=F,row.names=F,col.names=T,sep="\t")
