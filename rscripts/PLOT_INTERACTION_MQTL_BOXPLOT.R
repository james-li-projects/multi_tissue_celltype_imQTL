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
    write.table(data.frame(0,cov_df$Sample_Name),file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list2",quote=F,row.names=F,col.names=F,sep="\t")
    
    # iterating across all cell type combinations
    setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",dataset,"/imQTL/top_assoc"))
    tmp_combination_list <- list.files()[grepl(current_tissue,list.files()) & grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
    tmp_combination_list <- gsub("tensorQTL_imQTL_","",gsub("\\.cis_qtl_top_assoc.txt.gz","",tmp_combination_list))
    # removing eosinophils
    tmp_combination_list <- tmp_combination_list[!grepl("Eosino",tmp_combination_list)]
    for (tmp_combination in tmp_combination_list) {
      
      ###########################################
      ###########################################
      ###########################################
      # Extract unprocessed names
      unprocessed_tissue <- str_split(tmp_combination, pattern = "_")[[1]][1]
      unprocessed_celltype <- str_split(tmp_combination, pattern = "_")[[1]][2]
      # Create lookup vectors
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
      # Map to beautiful names
      processed_tissue <- tissue_map[[unprocessed_tissue]]
      processed_celltype <- celltype_map[[unprocessed_celltype]]
      processed_tmp_combination <- paste0(processed_tissue," - ",processed_celltype)
      ###########################################
      ###########################################
      ###########################################
      
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
        
        
        
        # Create an output directory for the tissue cell type combination if it does not exist
        dir_path2 <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/interaction_boxplots/", tmp_combination)
        if (!dir.exists(dir_path2)) {
          dir.create(dir_path2, recursive = TRUE)
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
          window_variants <- pvar %>% filter(`#CHROM`==paste0("chr",num_chr)) %>% filter(POS>(num_pos-25000),POS<(num_pos+25000))
          write.table(window_variants %>% select(ID),file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list2",quote=F,row.names=F,col.names=F,sep="\t")
          
          # filering DNAm bed file 
          dnam_bed <- input_dnam_bed %>% filter(phenotype_id==target_cpg) %>% head(1)
          t_dnam_bed <- data.frame(t(dnam_bed))
          t_dnam_bed <- t_dnam_bed %>% slice(-c(1:4)) 
          cpg_df <- t_dnam_bed %>% rename(mval=t.dnam_bed.) %>% mutate(Sample_Name=rownames(t_dnam_bed)) %>% select(Sample_Name,everything())
          
          # importing variant data
          system(paste0("module load plink/2.0; plink2 -pfile /gpfs/data/pierce-lab/james.li/imQTL/data/",dataset,"/genetic_data/processed_genetic_data_chrprefix --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list2 --keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list2 --maf 0.05 --export Av --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window2"))
          traw <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window2.traw") %>% separate(SNP,into=c("T1","T2","T3","T4"),remove=F) %>% mutate(across(11:ncol(.), ~ ifelse(COUNTED != T4, 2 - ., .))) %>% select(-c(paste0(rep("T","4"),1:4)))
          
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
          
          
          
          
          #####################################
          # preparing DF to plot interaction boxplot
          interaction_boxplot_df <- reg_df %>% select(mval,CT,all_of(lead_variant_period)) %>% rename_with(~ "dosage", .cols = 3) %>% mutate(dosage=round(dosage,digits=0)) %>% filter(dosage%in%c(0,1,2))
          
          hom_rec_str<-paste0(  
            unlist(str_split(lead_variant_period,pattern="\\."))[3]
            ,"/",
            unlist(str_split(lead_variant_period,pattern="\\."))[3]
          )
          
          het_str<-paste0(  
            unlist(str_split(lead_variant_period,pattern="\\."))[3]
            ,"/",
            unlist(str_split(lead_variant_period,pattern="\\."))[4]
          )
          hom_dom_str<-paste0(  
            unlist(str_split(lead_variant_period,pattern="\\."))[4]
            ,"/",
            unlist(str_split(lead_variant_period,pattern="\\."))[4]
          )
          
          interaction_boxplot_df <- interaction_boxplot_df %>% mutate(Genotype=ifelse(dosage==0,hom_rec_str,ifelse(dosage==1,het_str,ifelse(dosage==2,hom_dom_str,NA)))) %>% filter(!is.na(Genotype)) 
          
          # creating facets of upper/lower 50% and all individuals
          interaction_boxplot_df <- interaction_boxplot_df %>% mutate(`Cell Type Proportion`=ifelse(CT>median(CT),paste0("Upper 50% (",processed_celltype,")"),paste0("Lower 50% (",processed_celltype,")")))
          copy_interaction_boxplot_df<-interaction_boxplot_df %>% mutate(`Cell Type Proportion`="All Individuals")
          interaction_boxplot_df <- rbind(interaction_boxplot_df,copy_interaction_boxplot_df)
          
          # create numeric version of Genotype for regression lines in plot
          interaction_boxplot_df <- interaction_boxplot_df %>%
            mutate(Genotype_numeric = as.numeric(as.factor(Genotype)))
          
          # Plot
          interact_boxplot_p <- ggplot(interaction_boxplot_df, aes(x = Genotype, y = mval)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.5, aes(color = Genotype)) +  # visualize actual points
            geom_smooth(
              aes(x = Genotype_numeric, y = mval), 
              method = "lm", se = FALSE, color = "darkgray", linewidth = 1
            ) +
            facet_wrap(~ `Cell Type Proportion`, nrow = 1) +
            theme_classic() +
            labs(
              x = "Genotype",
              y = "DNAm (INT M-value)",
              title = paste(lead_variant_raw, target_cpg, sep = " - ")
            ) +
            theme(
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 16),
              strip.text = element_text(size = 14),
              plot.title = element_text(size = 16, face = "bold"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
              legend.position = "none"
            )
          
          # Save the plot
          output_path <- paste0(dir_path2, "/", lead_variant_period, "_", target_cpg, ".png")
          ggsave(output_path, interact_boxplot_p, width = 11, height = 5)
          
        }
      }
    }
  }
}
