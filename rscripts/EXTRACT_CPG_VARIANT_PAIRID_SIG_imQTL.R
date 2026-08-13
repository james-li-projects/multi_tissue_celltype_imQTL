library(data.table)
library(dplyr)
library(tidyr)

# extracting cpg and variant ids
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")
#all_pair_id <- wide_parsed_imqtl %>% filter(Category=="Same") %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_")) %>% select(pair_id) %>% unique()
all_pair_id <- wide_parsed_imqtl %>% mutate(pair_id=paste(phenotype_id,variant_id,sep="_")) %>% select(pair_id) %>% unique()
write.table(data.frame(all_pair_id),file="/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cpg_variant_pairs/sig_pair_id.list",quote=F,row.names=F,col.names=T)
