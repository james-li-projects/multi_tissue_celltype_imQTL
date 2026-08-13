library(data.table)
library(dplyr)

# setting import directory
setwd("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/eligible_gtex_samples")

# importing tissue argument
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
print(paste("PROCESSING INPUT FILES FOR:",tissue))

# importing bed file
bed <- fread(paste0(tissue,".bed.gz"))

# importing covariates file
cov <- fread(paste0(tissue,".covariates.txt"))

# importing psam file
psam <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix.psam"))
print(dim(bed))
print(dim(cov))

# identifying samples with complete data for mqtl mapping
eligible_sample_list <- intersect(intersect(colnames(bed),colnames(cov)),psam$`#IID`)

# function to convert tissue string to output string
convert_tissue <- function(input_tissue) {
  tissue_map <- list(
    "BreastMammaryTissue" = "breast",
    "ColonTransverse" = "colon",
    "KidneyCortex" = "kidney",
    "Lung" = "lung",
    "MuscleSkeletal" = "muscle",
    "Ovary" = "ovary",
    "Prostate" = "prostate",
    "Testis" = "testis",
    "WholeBlood" = "wb"
  )
  output_tissue <- tissue_map[[input_tissue]]
  if (is.null(output_tissue)) {
    stop("Invalid tissue type provided.")
  }
  return(output_tissue)
}

# save lists of complete data samples"
output_tissue <- convert_tissue(tissue)
save(eligible_sample_list,file=paste0("eligible_sample_list_",output_tissue,".RData"))
