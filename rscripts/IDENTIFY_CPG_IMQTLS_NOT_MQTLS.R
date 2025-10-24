# Ensure required libraries are loaded
library(dplyr)

# loading lists of CpGs
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_cpg_df.RData")
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/imQTL_cpg_df.RData")

# Create the list of unique tissues in imQTL
unique_tissues <- unique(imQTL_cpg_df$Tissue)

# Initialize empty list to store results
imQTL_specific_phenotypes <- list()

# Loop through each tissue
for (tissue in unique_tissues) {
  # Get phenotype IDs for this tissue from imQTL
  imQTL_ids <- imQTL_cpg_df %>%
    filter(Tissue == tissue) %>%
    pull(phenotype_id) %>%
    unique()
  
  # Get phenotype IDs for this tissue from mQTL
  mQTL_ids <- mQTL_cpg_df %>%
    filter(Tissue == tissue) %>%
    pull(phenotype_id) %>%
    unique()
  
  # Identify imQTL-specific phenotype IDs
  specific_ids <- setdiff(imQTL_ids, mQTL_ids)
  
  # Add to list
  imQTL_specific_phenotypes[[tissue]] <- specific_ids
}

# Convert list to data frame, skipping empty elements
imQTL_specific_df <- do.call(rbind, lapply(names(imQTL_specific_phenotypes), function(tissue) {
  ids <- imQTL_specific_phenotypes[[tissue]]
  if (length(ids) > 0) {
    data.frame(
      Tissue = tissue,
      phenotype_id = ids,
      stringsAsFactors = FALSE
    )
  }
}))

# View the result
head(imQTL_specific_df)

# writing out table
write.table(imQTL_specific_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/cpgs_imqtls_not_mqtls.tsv",quote=F,row.names=F,col.names=T,sep="\t")