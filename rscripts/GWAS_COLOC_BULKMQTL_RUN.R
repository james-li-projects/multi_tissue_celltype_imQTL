library(data.table)
library(coloc)

################################
# IMPORTING ARGUMENTS FOR GWAS #
################################
args <- commandArgs(trailingOnly = TRUE)
tmp_GWAS_file <- args[1]

#################################
# IDENTIFYING CASE-CONTROL GWAS #
#################################
GWAS_dir <- "/gpfs/data/pierce-lab/james.li/GWAS_Trait_87"

GWAS_list <- list.files(
  GWAS_dir,
  pattern = "\\.txt\\.gz$"
)

cc_traits <- unique(c(
  GWAS_list[grepl("self_report", GWAS_list)],
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
))

gwas_basename <- basename(tmp_GWAS_file)
gwas_type <- if (gwas_basename %chin% cc_traits) "cc" else "quant"

#################################
# IMPORTING AND PROCESSING GWAS #
#################################

# Import only columns actually needed
tmp_GWAS <- fread(
  tmp_GWAS_file,
  select = c(
    "chromosome",
    "position",
    "non_effect_allele",
    "effect_allele",
    "effect_size",
    "standard_error",
    "frequency",
    "sample_size"
  ),
  showProgress = TRUE
)

# Modify by reference rather than repeatedly copying the table
tmp_GWAS[
  ,
  `:=`(
    variant_id  = paste(
      chromosome,
      position,
      non_effect_allele,
      effect_allele,
      sep = ":"
    ),
    beta        = effect_size,
    varbeta     = standard_error^2,
    MAF         = pmin(frequency, 1 - frequency),
    N           = sample_size
  )
]

tmp_GWAS <- tmp_GWAS[
  complete.cases(variant_id, beta, varbeta, MAF, N) &
    MAF > 0 &
    MAF < 1 &
    varbeta > 0,
  .(variant_id, beta, varbeta, MAF, N)
]

# Ensure one row per variant
tmp_GWAS <- unique(tmp_GWAS, by = "variant_id")

# Critical speed improvement: index GWAS by variant ID
setkey(tmp_GWAS, variant_id)

# Calculate these once rather than for every mQTL
gwas_N <- median(tmp_GWAS$N, na.rm = TRUE)

rm_cols <- c(
  "chromosome",
  "position",
  "non_effect_allele",
  "effect_allele",
  "effect_size",
  "standard_error",
  "frequency",
  "sample_size"
)

gc()

#######################################
# IMPORTING LIST OF EXTRACTED WINDOWS #
#######################################
mqtl_dir <- "/scratch/jll1/imQTL/coloc_windows/extracted_mqtl_windows"

mqtl_file_list <- list.files(
  mqtl_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

total_num_mqtl <- length(mqtl_file_list)

##########################################
# INITIALIZING LIST TO STORE THE RESULTS #
##########################################

# Do not repeatedly rbind a growing data.frame
coloc_results_list <- vector("list", total_num_mqtl)
result_index <- 0L

#################################
# ITERATING THROUGH ALL IMQTLS  #
#################################
for (i in seq_along(mqtl_file_list)) {
  
  # Printing every iteration can be surprisingly expensive
  if (i == 1L || i %% 100L == 0L || i == total_num_mqtl) {
    message(
      "Performing coloc on mQTL ",
      i,
      " out of ",
      total_num_mqtl
    )
  }
  
  tmp_mqtl_file <- mqtl_file_list[i]
  tmp_mqtl_basename <- basename(tmp_mqtl_file)
  
  ####################################
  # PARSING INFORMATION FROM FILENAME #
  ####################################
  tmp_cleaned <- sub("\\.tsv$", "", tmp_mqtl_basename)
  parts <- strsplit(tmp_cleaned, "_", fixed = TRUE)[[1]]
  
  tmp_mqtl_combination <- parts[1]
  tmp_mqtl_cpg_id <- parts[2]
  tmp_mqtl_variant_id <- paste(parts[3:6], collapse = ":")
  
  ##################################
  # IMPORTING AND PROCESSING IMQTL #
  ##################################
  
  # Read only required columns
  tmp_mqtl <- fread(
    tmp_mqtl_file,
    select = c(
      "variant_id",
      "af",
      "ma_count",
      "slope",
      "slope_se"
    ),
    showProgress = FALSE
  )
  
  tmp_mqtl[
    ,
    `:=`(
      beta    = slope,
      varbeta = slope_se^2,
      MAF     = pmin(af, 1 - af),
      N       = ma_count / (2 * af)
    )
  ]
  
  tmp_mqtl <- tmp_mqtl[
    complete.cases(variant_id, beta, varbeta, MAF, N) &
      MAF > 0 &
      MAF < 1 &
      varbeta > 0,
    .(variant_id, beta, varbeta, MAF, N)
  ]
  
  tmp_mqtl <- unique(tmp_mqtl, by = "variant_id")
  
  # Check the lead variant in the mQTL window without scanning the GWAS
  if (!(tmp_mqtl_variant_id %chin% tmp_mqtl$variant_id)) {
    next
  }
  
  ######################################
  # FAST KEYED JOIN WITH GWAS DATASET  #
  ######################################
  
  # This replaces:
  # intersect()
  # two %in% filters
  # two separate dataset constructions
  #
  # nomatch = 0L retains only shared variants
  common <- tmp_GWAS[
    tmp_mqtl,
    on = "variant_id",
    nomatch = 0L,
    .(
      snp = variant_id,
      gwas_beta = x.beta,
      gwas_varbeta = x.varbeta,
      gwas_MAF = x.MAF,
      mqtl_beta = i.beta,
      mqtl_varbeta = i.varbeta,
      mqtl_MAF = i.MAF,
      mqtl_N = i.N
    )
  ]
  
  # Need >50 shared variants and the lead variant in the overlap
  if (
    nrow(common) <= 50L ||
    !(tmp_mqtl_variant_id %chin% common$snp)
  ) {
    next
  }
  
  #################################
  # CONSTRUCTING COLOC DATASETS   #
  #################################
  coloc_dataset1 <- list(
    snp     = common$snp,
    beta    = common$gwas_beta,
    varbeta = common$gwas_varbeta,
    MAF     = common$gwas_MAF,
    N       = gwas_N,
    type    = gwas_type
  )
  
  coloc_dataset2 <- list(
    snp     = common$snp,
    beta    = common$mqtl_beta,
    varbeta = common$mqtl_varbeta,
    MAF     = common$mqtl_MAF,
    N       = median(common$mqtl_N, na.rm = TRUE),
    type    = "quant"
  )
  
  ####################
  # PERFORMING COLOC #
  ####################
  coloc_output <- coloc.abf(
    dataset1 = coloc_dataset1,
    dataset2 = coloc_dataset2
  )
  
  coloc_summary <- coloc_output$summary
  
  result_index <- result_index + 1L
  
  coloc_results_list[[result_index]] <- data.table(
    combination = tmp_mqtl_combination,
    cpg_id      = tmp_mqtl_cpg_id,
    variant_id  = tmp_mqtl_variant_id,
    trait       = tmp_GWAS_file,
    PP.H0.abf   = unname(coloc_summary["PP.H0.abf"]),
    PP.H1.abf   = unname(coloc_summary["PP.H1.abf"]),
    PP.H2.abf   = unname(coloc_summary["PP.H2.abf"]),
    PP.H3.abf   = unname(coloc_summary["PP.H3.abf"]),
    PP.H4.abf   = unname(coloc_summary["PP.H4.abf"])
  )
}

################################
# COMBINING ALL COLOC RESULTS  #
################################
if (result_index > 0L) {
  coloc_results <- rbindlist(
    coloc_results_list[seq_len(result_index)],
    use.names = TRUE
  )
} else {
  coloc_results <- data.table(
    combination = character(),
    cpg_id      = character(),
    variant_id  = character(),
    trait       = character(),
    PP.H0.abf   = numeric(),
    PP.H1.abf   = numeric(),
    PP.H2.abf   = numeric(),
    PP.H3.abf   = numeric(),
    PP.H4.abf   = numeric()
  )
}

############################
# OUTPUTTING COLOC RESULTS #
############################
output_dir <- paste0(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/",
  "GWAS_coloc_bulkmqtl/coloc_output"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  output_dir,
  paste0(sub("\\.txt\\.gz$", "", gwas_basename), ".tsv")
)

fwrite(
  coloc_results,
  file = output_file,
  sep = "\t",
  quote = FALSE
)

message("Completed. Results written to: ", output_file)