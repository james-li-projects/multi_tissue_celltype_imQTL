library(dplyr)
library(data.table)
library(coloc)

################################
# IMPORTING ARGUMENTS FOR GWAS #
################################
# importing file arguments for GWAS
args <- commandArgs(trailingOnly = TRUE)
tmp_GWAS_file <- args[1]

#################################
# IMPORTING AND PROCESSING GWAS #
#################################
# obtaining list of all potential GWAS traits
GWAS_list=list.files("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87")[grepl(".txt.gz",list.files("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87"))]

# identifying which GWAS traits are case-control studies
cc_traits <- c()
cc_traits <- c(cc_traits,GWAS_list[grepl("self_report",GWAS_list)])
cc_traits <- c(
  cc_traits,
  "imputed_BCAC_ER_negative_BreastCancer_EUR.txt.gz",
  "imputed_BCAC_ER_positive_BreastCancer_EUR.txt.gz",
  "imputed_BCAC_Overall_BreastCancer_EUR.txt.gz",
  "imputed_CNCR_Insomnia_all.txt.gz",
  "imputed_EAGLE_Eczema.txt.gz",
  "imputed_IBD.EUR.Crohns_Disease.txt.gz",
  "imputed_IBD.EUR.Inflammatory_Bowel_Disease.txt.gz",
  "imputed_IBD.EUR.Ulcerative_Colitis.txt.gz",
  "imputed_IGAP_Alzheimer.txt.gz",
  "imputed_IMMUNOBASE_Systemic_lupus_erythematosus_hg19.txt.gz",
  "imputed_Jones_et_al_2016_Chronotype.txt.gz",
  "imputed_PGC_ADHD_EUR_2017.txt.gz",
  "imputed_RA_OKADA_TRANS_ETHNIC.txt.gz",
  "imputed_SSGAC_Depressive_Symptoms.txt.gz",
  "imputed_UKB_1180_Morning_or_evening_person_chronotype.txt.gz",
  "imputed_UKB_1200_Sleeplessness_or_insomnia.txt.gz",
  "imputed_UKB_2395_2_Hair_or_balding_pattern_Pattern_2.txt.gz",
  "imputed_UKB_2395_3_Hair_or_balding_pattern_Pattern_3.txt.gz",
  "imputed_UKB_2395_4_Hair_or_balding_pattern_Pattern_4.txt.gz",
  "imputed_UKB_6150_1_Vascular_or_heart_problems_diagnosed_by_doctor_Heart_attack.txt.gz",
  "imputed_UKB_6152_5_diagnosed_by_doctor_Blood_clot_in_the_leg_DVT.txt.gz",
  "imputed_UKB_6152_7_diagnosed_by_doctor_Blood_clot_in_the_lung.txt.gz",
  "imputed_UKB_6152_8_diagnosed_by_doctor_Asthma.txt.gz",
  "imputed_UKB_6152_9_diagnosed_by_doctor_Hayfever_allergic_rhinitis_or_eczema.txt.gz",
  "imputed_UKB_G40_Diagnoses_main_ICD10_G40_Epilepsy.txt.gz",
  "imputed_UKB_G43_Diagnoses_main_ICD10_G43_Migraine.txt.gz",
  "imputed_pgc.scz2.txt.gz",
  "cancer_colon.txt.gz",
  "cancer_leukemia.txt.gz",
  "cancer_lung.txt.gz",
  "cancer_ovary.txt.gz",
  "cancer_prostate.txt.gz"
)
cc_traits <- unique(cc_traits)

# importing and parsing GWAS DF
tmp_GWAS <- fread(tmp_GWAS_file) %>% mutate(trait=tmp_GWAS_file)
tmp_GWAS <- tmp_GWAS %>% mutate(variant_id = paste(chromosome,position,non_effect_allele,effect_allele,sep=":")) %>% mutate(GWAS_varbeta = standard_error^2,GWAS_MAF=ifelse(frequency > 0.5,1-frequency,frequency),GWAS_type=ifelse(trait %in% cc_traits,"cc","quant")) %>% rename(GWAS_beta=effect_size,GWAS_N=sample_size) %>% select(variant_id,GWAS_beta,GWAS_varbeta,GWAS_MAF,GWAS_N,trait,GWAS_type) %>% filter(!is.na(GWAS_beta))
tmp_GWAS <- na.omit(tmp_GWAS)

##########################################
# INITIALIZING DF TO STORE COLOC RESULTS #
##########################################
# initializing df to store coloc results
coloc_results <-data.frame()

#################################
# IMPORTING ARGUMENTS FOR IMQTL #
#################################
# set working directory to imQTL extracted data files
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_data")
imqtl_file_list <- list.files()
total_num_imqtl <- length(imqtl_file_list)

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
  tmp_imqtl_cpg_id <- parts[3]
  tmp_imqtl_variant_id <- paste(parts[4:7], collapse = ":")
  
  ##################################
  # IMPORTING AND PROCESSING IMQTL #
  ##################################
  # importing imqtl DF
  tmp_imqtl <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_data/",tmp_imqtl_file))
  tmp_imqtl <- tmp_imqtl %>% mutate(imqtl_N=ma_count/(2*af),imqtl_MAF=ifelse(af>0.5,1-af,af),imqtl_varbeta=b_gi_se^2,imqtl_beta=b_gi,imqtl_type="quant") %>% select(variant_id,imqtl_beta,imqtl_varbeta,imqtl_MAF,imqtl_N,phenotype_id,imqtl_type)
  tmp_imqtl <- na.omit(tmp_imqtl)
  
  ####################
  # PERFORMING COLOC #
  ####################
  # identifying shared variants between both datasets
  tmp_common_variant_list <- intersect(tmp_GWAS$variant_id,tmp_imqtl$variant_id)
  
  # if more than 50 variants shared between both datasets continue processing for running coloc
  if (length(tmp_common_variant_list) > 50) {
    
    # creating coloc datasets
    dataset1 <- tmp_GWAS %>% filter(variant_id %in% tmp_common_variant_list) %>% rename(snp=variant_id) %>% select(snp,GWAS_beta,GWAS_varbeta,GWAS_MAF,GWAS_N,GWAS_type,trait)
    colnames(dataset1) <- gsub("GWAS_","",colnames(dataset1))
    dataset2 <- tmp_imqtl %>% filter(variant_id %in% tmp_common_variant_list) %>% rename(snp=variant_id) %>% select(snp,imqtl_beta,imqtl_varbeta,imqtl_MAF,imqtl_N,imqtl_type,phenotype_id)
    colnames(dataset2) <- gsub("imqtl_","",colnames(dataset2))
    
    # Convert dataset1 to list
    coloc_dataset1 <- list(
      snp     = dataset1$snp,
      beta    = dataset1$beta,
      varbeta = dataset1$varbeta,
      MAF     = dataset1$MAF,
      N       = median(unique(dataset1$N)),  # Ensure it's a single number
      type    = unique(dataset1$type)  # Should be "quant"
    )
    
    # Convert dataset2 to list
    coloc_dataset2 <- list(
      snp     = dataset2$snp,
      beta    = dataset2$beta,
      varbeta = dataset2$varbeta,
      MAF     = dataset2$MAF,
      N       = median(unique(dataset2$N)),
      type    = unique(dataset2$type)
    )
    
    # running coloc and extracting PP.H4
    coloc_output <- coloc.abf(
      coloc_dataset1,
      coloc_dataset2
    )
    PP.H0.abf=coloc_output$summary["PP.H0.abf"]
    PP.H1.abf=coloc_output$summary["PP.H1.abf"]
    PP.H2.abf=coloc_output$summary["PP.H2.abf"]
    PP.H3.abf=coloc_output$summary["PP.H3.abf"]
    PP.H4.abf=coloc_output$summary["PP.H4.abf"]
    
    # storing results
    tmp_coloc_results <- data.frame(
      combination=tmp_imqtl_combination,
      cpg_id=tmp_imqtl_cpg_id,
      variant_id=tmp_imqtl_variant_id,
      trait=tmp_GWAS_file,
      PP.H0.abf=PP.H0.abf,
      PP.H1.abf=PP.H1.abf,
      PP.H2.abf=PP.H2.abf,
      PP.H3.abf=PP.H3.abf,
      PP.H4.abf=PP.H4.abf
    )
    coloc_results <- rbind(
      coloc_results,
      tmp_coloc_results
    ) 
  }
}

############################
# OUTPUTTING COLOC RESULTS #
############################
# writing out coloc results as a tsv
write.table(coloc_results,file=paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/coloc_output/",gsub("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/","",gsub(".txt.gz","",tmp_GWAS_file)),".tsv"),quote=F,row.names=F,col.names=T,sep="\t")
