library(data.table)
library(dplyr)
extracted_dictionary_df<-fread("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/GWAS_Trait_87_Dictionary/PMID_33499903_GWAS_METADATA.txt") %>% mutate(V2=Phenotype,V3=paste0("imputed_",Tag,".txt.gz")) %>% select(V2,V3)

# modifying breast cancer columns
extracted_dictionary_df <- extracted_dictionary_df %>%
  mutate(
    V2 = case_when(
      V2 == "Estrogen-receptor-positive breast cancer in Europeans, imputed genotype" ~
        "Breast Cancer - Estrogen-receptor positive - Only Europeans",
      
      V2 == "Estrogen-receptor-negative breast cancer in Europeans, imputed genotype" ~
        "Breast Cancer - Estrogen-receptor negative - Only Europeans",
      
      V2 == "Overall breast cancer in Europeans, imputed genotype" ~
        "Breast Cancer (Overall) - Only Europeans",
      
      V2 == "Coronary Artery Disease (additive model)" ~
        "Coronary Artery Disease",
      
      V2 == "Attention deficit hyperactivity disorder - EUR" ~
        "Attention deficit - Hiperactivity disorder - Only Europeans",
      TRUE ~ V2
    )
  )

output_file <- paste0(
  "/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/GWAS_Trait_87_Dictionary/GWASes.87.txt"
)

######################################
# Write headerless tab-delimited file
######################################

fwrite(
  extracted_dictionary_df,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)

cat("\nCreated:", output_file, "\n")
