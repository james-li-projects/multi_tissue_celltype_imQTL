library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#######################################
# extracting numbers of imQTLs
Dataset <- c()
combination <- c()

imqtl_cpg_variant_df <- data.frame()

for (Dataset in c("GTEx","HEALS")) {
  setwd(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",Dataset,"/imQTL/top_assoc"))
  file_list <- list.files()[grepl(".cis_qtl_top_assoc.txt.gz",list.files())]
  for (combination in file_list) {
    print(combination)
    tmp_df <- fread(combination)
    # filtering for Bonferroni-adjusted p-value<0.05
    tmp_df <- tmp_df %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni"))
    tmp_df <- tmp_df %>% filter(pval_adj_bonf < 0.05)

    # assembling a DF containing all the cpgs
    tmp_imqtl_cpg_variant_df <- tmp_df %>% select(phenotype_id,variant_id) %>% mutate(combination=combination) %>% mutate(Dataset=Dataset)
    imqtl_cpg_variant_df <- rbind(imqtl_cpg_variant_df,tmp_imqtl_cpg_variant_df)
  }
}

# parsing out chromosome and combination
parsed_imqtl_cpg_variant_df <- imqtl_cpg_variant_df %>% separate(variant_id,remove=F,into=c("chr","pos","a2","a1"),sep="\\:") %>% mutate(chunk_num=as.integer(gsub("chr","",chr))) %>% mutate(parsed_combination=gsub("tensorQTL_imQTL_","",gsub("\\.cis_qtl_top_assoc.txt.gz","",combination))) %>% mutate(input_parquet_filename = paste0(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/",
  Dataset,
  "/imQTL/full_assoc/tensorQTL_imQTL_",
  parsed_combination,
  ".chunk",
  chunk_num,
  ".cis_qtl_pairs.chr",
  chunk_num,
  ".parquet"
)) %>% select(-a2,-a1,-combination,-Dataset)
# filtering out eosinophils
parsed_imqtl_cpg_variant_df <- parsed_imqtl_cpg_variant_df %>% filter(parsed_combination!="wb_Eosino")

##########################################
# dividing up this DF into 50 split DFs with approximately equal numbers of rows
output_dir_prefix <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_list/query_"
# Number of splits
n_splits <- 50
# Split into 50 approximately equal-sized parts
split_list <- split(parsed_imqtl_cpg_variant_df, cut(seq_len(nrow(parsed_imqtl_cpg_variant_df)), n_splits, labels = FALSE))
# Write each split to a separate TSV file
for (i in seq_along(split_list)) {
  output_path <- paste0(output_dir_prefix, i, ".tsv")
  write.table(split_list[[i]], file = output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
