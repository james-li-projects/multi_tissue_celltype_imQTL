library(data.table)
library(dplyr)
library(tidyr)

load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")

# Define lookup maps
tissue_map <- c(
  "breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
  "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary"
)

celltype_map <- c(
  "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
  "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
  "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
  "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
  "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
  "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
  "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
)

# Apply the mappings
wide_parsed_imqtl$tissue <- tissue_map[wide_parsed_imqtl$tissue]
wide_parsed_imqtl$celltype <- celltype_map[wide_parsed_imqtl$celltype]
write.table(wide_parsed_imqtl,file="/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl_CLEAN_CT_TISSUE_NAMES.tsv",quote=F,row.names=F,col.names=T,sep="\t")


# writing out cell type specific sumstats
ct_spec_imqtl_df <-fread("/gpfs/data/pierce-lab/james.li/imQTL/output/CANDIDATE_CT_SPECIFIC/tables/candidate_ct_specific_imqtl_DF.tsv") %>% mutate(tissue = recode(tissue,"colon" = "Colon","wb" = "Whole Blood")) %>% select(-dataset) %>% rename(variant_id=lead_variant_raw,phenotype_id=target_cpg)
ct_sumstats <- inner_join(wide_parsed_imqtl,ct_spec_imqtl_df)
write.table(ct_sumstats,file="/gpfs/data/pierce-lab/james.li/imQTL/output/CANDIDATE_CT_SPECIFIC/tables/ct_spec_sumstats.tsv",row.names=F,col.names=T,quote=F,sep="\t")
