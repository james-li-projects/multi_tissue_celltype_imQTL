library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# importing query file name
args <- commandArgs(trailingOnly = TRUE)
query_file <- args[1]

# importing query file
extract_list <- fread(query_file) %>% separate(parsed_combination,into=c("tissue","celltype"),sep="_",remove=F) %>%
  mutate(
    input_tsv_filename = input_parquet_filename %>%
      str_replace("tensorQTL_imQTL", "tensorQTL_mQTL") %>%
      str_replace(parsed_combination, tissue) %>%
      str_replace("imQTL/full_assoc", "mQTL/full_assoc_tsv") %>%
      str_replace("\\.parquet$", ".tsv") %>%
      str_replace(paste0("chunk",chunk_num,"."),"")
  ) %>% select(phenotype_id,variant_id,chr,pos,chunk_num,celltype,tissue,input_tsv_filename)

# parsing query file
setDT(extract_list)
extract_list[, c("chr_ex", "pos_ex", "ref_ex", "alt_ex") := tstrsplit(variant_id, ":", fixed = TRUE)]
extract_list[, pos_ex := as.integer(pos_ex)]

# Group extract_list by input file and process efficiently
input_file_groups <- split(extract_list, by = "input_tsv_filename")

# Loop over each group (i.e., each unique input file)
for (file_path in names(input_file_groups)) {
  # Read in the mQTL file just once
  tmp_mqtl_df <- fread(file_path)
  setDT(tmp_mqtl_df)
  
  # Parse variant_id
  tmp_mqtl_df[, c("chr", "pos", "ref", "alt") := tstrsplit(variant_id, ":", fixed = TRUE)]
  tmp_mqtl_df[, pos := as.integer(pos)]
  
  # Get the relevant subset of extract_list for this file
  sub_extract <- input_file_groups[[file_path]]
  
  # Loop over each row in the subset
  for (j in seq_len(nrow(sub_extract))) {
    row <- sub_extract[j]
    
    # Filter the mQTL data
    tmp_filtered_mqtl_df <- tmp_mqtl_df[
      phenotype_id == row[["phenotype_id"]] &
        chr == row[["chr_ex"]] &
        pos >= (row[["pos_ex"]] - 250000) &
        pos <= (row[["pos_ex"]] + 250000)
    ]
    
    # Define output path
    output_path <- paste0(
      "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_mqtl_data/",
      paste(row$tissue, row$celltype, row$phenotype_id, gsub(":", "_", row$variant_id), sep = "_"),
      ".tsv"
    )
    
    # Write filtered results
    fwrite(tmp_filtered_mqtl_df, file = output_path, sep = "\t", quote = FALSE)
  }
}
