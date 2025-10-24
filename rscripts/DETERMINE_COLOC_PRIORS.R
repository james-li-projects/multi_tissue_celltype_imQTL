library(data.table)
library(dplyr)
library(topr)

load("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/GWAS_DF.RData")

#############################
tmp_input <- GWAS_DF %>% mutate(CHROM=chromosome,REF=non_effect_allele,ALT=effect_allele,ID=panel_variant_id,BETA=effect_size,SE=standard_error,P=pvalue,POS=position) %>% select(CHROM,POS,ID,REF,ALT,BETA,SE,P) %>% mutate(P=ifelse(P<3e-303,3e-303,P))

lead_variants <- get_lead_snps(tmp_input,thresh = 5e-08,region_size = 1e+06)
num_lead_GWAS_variants <- nrow(lead_variants)
num_total_GWAS_variants <- 9000000
num_lead_GWAS_variants/num_total_GWAS_variants

#############################
mqtl_variant_bed<-fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed/all_variant_GTEx_lung.bed")



