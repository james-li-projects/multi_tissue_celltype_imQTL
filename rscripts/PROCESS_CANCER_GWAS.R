library(data.table)
library(dplyr)
library(tidyr)

# importing allele frequencies from GTEx and HEALS
afreq_gtex <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix_FREQ.afreq")
afreq_heals <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/genetic_data/processed_genetic_data_chrprefix_FREQ.afreq")
afreq_combined <- rbind(
  afreq_gtex,
  afreq_heals %>% filter(!(ID%in%afreq_gtex$ID))
)
afreq_combined <- afreq_combined %>% mutate(panel_variant_id=paste(gsub("\\:","_",ID),"b38",sep="_")) %>% mutate(MAF=ifelse(ALT_FREQS>0.5,1-ALT_FREQS,ALT_FREQS)) %>% select(panel_variant_id,MAF)

################################
# Processing colon cancer GWAS # 
################################
# importing sumstats
input_sumstats <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GWAS_CATALOG/cancer/colon/GCST90129505_buildGRCh37.tsv")

# computing effective sample size
cases1 <- 21731
controls1 <- 47444
cases2 <- 78473
controls2 <- 107143
cases=cases1+cases2
controls=controls1+controls2

# Proportion of cases
s <- cases / (cases + controls)
# Effective sample size
neff <- 4 / (1 / cases + 1 / controls)
# Effective number of cases
neff_cases <- round(neff*s)

# parsing sumstats
#### parsed_input_sumstats <- input_sumstats %>% head()
parsed_input_sumstats <- input_sumstats %>% mutate(
  chromosome=paste0("chr",chromosome),
  position=as.integer(base_pair_location),
  initial_variant_id=Variant_id,
  effect_allele=effect_allele,
  non_effect_allele=other_allele,
  effect_size=beta,
  standard_error=standard_error,
  pvalue=p_value,
  sample_size=round(neff),
  n_cases=round(neff_cases)
)
parsed_input_sumstats <- parsed_input_sumstats %>% select(chromosome,position,initial_variant_id,effect_allele,non_effect_allele,effect_size,standard_error,pvalue,sample_size,n_cases)

tmp_liftOver <- parsed_input_sumstats %>% mutate(position2=as.integer(position+1)) %>% select(chromosome,position,position2,initial_variant_id)

# writing out and performing liftOver
setwd("/gpfs/data/pierce-lab/james.li/imQTL/tmp/liftOver")
write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
# reading in liftOver mapped results
liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
colnames(liftOver_mapped) <- c("chromosome","position_hg38","position_hg38_2","initial_variant_id")
liftOver_mapped <- liftOver_mapped %>% select(-position_hg38_2)

# joining liftOver coordinates - orientation #1
output_sumstats1 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,non_effect_allele,effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats1 <- output_sumstats1 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# joining liftOver coordinates - orientation #2
output_sumstats2 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,effect_allele,non_effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id,old_effect=effect_allele,old_non=non_effect_allele) %>% rename(effect_allele=old_non,non_effect_allele=old_effect) %>% mutate(effect_size = -effect_size) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats2 <- output_sumstats2 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# combining all variants
output_sumstats <- rbind(output_sumstats1,output_sumstats2)
                                                
# performing final parsing of sumstats to match the 87 trait GWAS sumstats
output_sumstats <- output_sumstats %>% mutate(
  variant_id=panel_variant_id, 
  current_build="hg38",
  frequency=MAF,
  imputation_status="imputed"
  ) %>% select(
    variant_id,panel_variant_id,chromosome,position,effect_allele,non_effect_allele,current_build,frequency,sample_size,zscore,pvalue,effect_size,standard_error,imputation_status,n_cases
    )

# removing records with na values
output_sumstats <- na.omit(output_sumstats)

# save only the first instance of any duplicated initial_variant_id
output_sumstats <- output_sumstats[!duplicated(panel_variant_id)]

# writing out final output summary statistics
write.table(
  output_sumstats,
  file = gzfile("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/cancer_colon.txt.gz"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

###############################
# Processing lung cancer GWAS # 
###############################
# importing sumstats
input_sumstats <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GWAS_CATALOG/cancer/lung/GCST90134661_buildGRCh37.tsv")

# computing effective sample size
cases1 <- 26683 
controls1 <- 25278 
cases2 <- 7062 
controls2 <- 5372 
cases3 <- 1987 
controls3 <- 3774 
cases=cases1+cases2+cases3
controls=controls1+controls2+controls3

# Proportion of cases
s <- cases / (cases + controls)
# Effective sample size
neff <- 4 / (1 / cases + 1 / controls)
# Effective number of cases
neff_cases <- round(neff*s)

# parsing sumstats
#### parsed_input_sumstats <- input_sumstats %>% head()
parsed_input_sumstats <- input_sumstats %>% mutate(
  chromosome=paste0("chr",chromosome),
  position=as.integer(base_pair),
  initial_variant_id=variant_id,
  effect_allele=effect_allele,
  non_effect_allele=other_allele,
  effect_size=log(odds_ratio),
  standard_error=standard_error,
  pvalue=p_value,
  sample_size=round(neff),
  n_cases=round(neff_cases)
)
parsed_input_sumstats <- parsed_input_sumstats %>% select(chromosome,position,initial_variant_id,effect_allele,non_effect_allele,effect_size,standard_error,pvalue,sample_size,n_cases)
parsed_input_sumstats <- na.omit(parsed_input_sumstats)
# save only the first instance of any duplicated initial_variant_id
parsed_input_sumstats <- parsed_input_sumstats[!duplicated(initial_variant_id)]

# preparing liftOver input table
tmp_liftOver <- parsed_input_sumstats %>% mutate(position2=as.integer(position+1)) %>% select(chromosome,position,position2,initial_variant_id)

# writing out and performing liftOver
setwd("/gpfs/data/pierce-lab/james.li/imQTL/tmp/liftOver")
write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
# reading in liftOver mapped results
liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
colnames(liftOver_mapped) <- c("chromosome","position_hg38","position_hg38_2","initial_variant_id")
liftOver_mapped <- liftOver_mapped %>% select(-position_hg38_2)

# joining liftOver coordinates - orientation #1
output_sumstats1 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,non_effect_allele,effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats1 <- output_sumstats1 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# joining liftOver coordinates - orientation #2
output_sumstats2 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,effect_allele,non_effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id,old_effect=effect_allele,old_non=non_effect_allele) %>% rename(effect_allele=old_non,non_effect_allele=old_effect) %>% mutate(effect_size = -effect_size) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats2 <- output_sumstats2 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# combining all variants
output_sumstats <- rbind(output_sumstats1,output_sumstats2)

# performing final parsing of sumstats to match the 87 trait GWAS sumstats
output_sumstats <- output_sumstats %>% mutate(
  variant_id=panel_variant_id, 
  current_build="hg38",
  frequency=MAF,
  imputation_status="imputed"
) %>% select(
  variant_id,panel_variant_id,chromosome,position,effect_allele,non_effect_allele,current_build,frequency,sample_size,zscore,pvalue,effect_size,standard_error,imputation_status,n_cases
)

# removing records with na values
output_sumstats <- na.omit(output_sumstats)

# save only the first instance of any duplicated initial_variant_id
output_sumstats <- output_sumstats[!duplicated(panel_variant_id)]

# writing out final output summary statistics
write.table(
  output_sumstats,
  file = gzfile("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/cancer_lung.txt.gz"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

############################
# Processing leukemia GWAS # 
############################
# importing sumstats
input_sumstats <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GWAS_CATALOG/cancer/leukemia/GCST90479830.tsv.gz")

# computing effective sample size
cases1 <- 360  
controls1 <- 121531  
cases2 <- 2845  
controls2 <- 448087  
cases=cases1+cases2
controls=controls1+controls2

# Proportion of cases
s <- cases / (cases + controls)
# Effective sample size
neff <- 4 / (1 / cases + 1 / controls)
# Effective number of cases
neff_cases <- round(neff*s)

# parsing sumstats
#### parsed_input_sumstats <- input_sumstats %>% head()
input_sumstats <- input_sumstats %>%
  mutate(
    odds_ratio = as.numeric(odds_ratio),
    ci_upper = as.numeric(ci_upper),
    ci_lower = as.numeric(ci_lower),
    effect_size = log(odds_ratio),
    standard_error = (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
  )
parsed_input_sumstats <- input_sumstats %>% mutate(
  chromosome=paste0("chr",chromosome),
  position=as.integer(base_pair_location),
  initial_variant_id=rsid,
  effect_allele=effect_allele,
  non_effect_allele=other_allele,
  pvalue=as.numeric(p_value),
  sample_size=round(neff),
  n_cases=round(neff_cases)
)
# selecting relevant columns
parsed_input_sumstats <- parsed_input_sumstats %>% select(chromosome,position,initial_variant_id,effect_allele,non_effect_allele,effect_size,standard_error,pvalue,sample_size,n_cases)
# removing records with na values
parsed_input_sumstats <- na.omit(parsed_input_sumstats)
# save only the first instance of any duplicated initial_variant_id
parsed_input_sumstats <- parsed_input_sumstats[!duplicated(initial_variant_id)]

# further parsing the variants
output_sumstats <- parsed_input_sumstats %>% mutate(panel_variant_id=paste(chromosome,position,effect_allele,non_effect_allele,"b38",sep="_")) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id,old_effect=effect_allele,old_non=non_effect_allele) %>% rename(effect_allele=old_non,non_effect_allele=old_effect) %>% mutate(effect_size = -effect_size) %>% mutate(zscore=effect_size/standard_error)

# joining allele frequencies
output_sumstats <- output_sumstats %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# performing final parsing of sumstats to match the 87 trait GWAS sumstats
output_sumstats <- output_sumstats %>% mutate(
  variant_id=panel_variant_id, 
  current_build="hg38",
  frequency=MAF,
  imputation_status="imputed"
) %>% select(
  variant_id,panel_variant_id,chromosome,position,effect_allele,non_effect_allele,current_build,frequency,sample_size,zscore,pvalue,effect_size,standard_error,imputation_status,n_cases
)

# removing records with na values
output_sumstats <- na.omit(output_sumstats)

# save only the first instance of any duplicated initial_variant_id
output_sumstats <- output_sumstats[!duplicated(panel_variant_id)]

# writing out final output summary statistics
write.table(
  output_sumstats,
  file = gzfile("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/cancer_leukemia.txt.gz"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

################################
# Processing ovary cancer GWAS # 
################################
# importing sumstats
input_sumstats <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GWAS_CATALOG/cancer/ovary/GCST90244168_buildGRCh38.tsv.gz")

# computing effective sample size
cases1 <- 15588  
controls1 <- 105724  
cases2 <- 0 
controls2 <- 0 
cases3 <- 0 
controls3 <- 0 
cases=cases1+cases2+cases3
controls=controls1+controls2+controls3

# Proportion of cases
s <- cases / (cases + controls)
# Effective sample size
neff <- 4 / (1 / cases + 1 / controls)
# Effective number of cases
neff_cases <- round(neff*s)

# parsing sumstats
#### parsed_input_sumstats <- input_sumstats %>% head()
parsed_input_sumstats <- input_sumstats %>% mutate(
  chromosome=paste0("chr",chromosome),
  position=as.integer(base_pair_location),
  initial_variant_id=SNP,
  effect_allele=effect_allele,
  non_effect_allele=other_allele,
  effect_size=beta,
  standard_error=standard_error,
  pvalue=p_value,
  sample_size=round(neff),
  n_cases=round(neff_cases)
)
parsed_input_sumstats <- parsed_input_sumstats %>% select(chromosome,position,initial_variant_id,effect_allele,non_effect_allele,effect_size,standard_error,pvalue,sample_size,n_cases)
parsed_input_sumstats <- na.omit(parsed_input_sumstats)
# save only the first instance of any duplicated initial_variant_id
parsed_input_sumstats <- parsed_input_sumstats[!duplicated(initial_variant_id)]

# preparing liftOver input table
tmp_liftOver <- parsed_input_sumstats %>% mutate(position2=as.integer(position+1)) %>% select(chromosome,position,position2,initial_variant_id)

# writing out and performing liftOver
setwd("/gpfs/data/pierce-lab/james.li/imQTL/tmp/liftOver")
write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
# reading in liftOver mapped results
liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
colnames(liftOver_mapped) <- c("chromosome","position_hg38","position_hg38_2","initial_variant_id")
liftOver_mapped <- liftOver_mapped %>% select(-position_hg38_2)

# joining liftOver coordinates - orientation #1
output_sumstats1 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,non_effect_allele,effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats1 <- output_sumstats1 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# joining liftOver coordinates - orientation #2
output_sumstats2 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,effect_allele,non_effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id,old_effect=effect_allele,old_non=non_effect_allele) %>% rename(effect_allele=old_non,non_effect_allele=old_effect) %>% mutate(effect_size = -effect_size) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats2 <- output_sumstats2 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# combining all variants
output_sumstats <- rbind(output_sumstats1,output_sumstats2)

# performing final parsing of sumstats to match the 87 trait GWAS sumstats
output_sumstats <- output_sumstats %>% mutate(
  variant_id=panel_variant_id, 
  current_build="hg38",
  frequency=MAF,
  imputation_status="imputed"
) %>% select(
  variant_id,panel_variant_id,chromosome,position,effect_allele,non_effect_allele,current_build,frequency,sample_size,zscore,pvalue,effect_size,standard_error,imputation_status,n_cases
)

# removing records with na values
output_sumstats <- na.omit(output_sumstats)

# save only the first instance of any duplicated initial_variant_id
output_sumstats <- output_sumstats[!duplicated(panel_variant_id)]

# writing out final output summary statistics
write.table(
  output_sumstats,
  file = gzfile("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/cancer_ovary.txt.gz"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

###################################
# Processing prostate cancer GWAS # 
###################################
# importing sumstats
input_sumstats <- fread("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GWAS_CATALOG/cancer/prostate/GCST90274713.tsv")

# computing effective sample size
cases1 <- 122188 
controls1 <- 604640
cases2 <- 19391 
controls2 <- 61608 
cases3 <- 10809 
controls3 <- 95790 
cases4 <- 3931 
controls4 <- 26405 
cases=cases1+cases2+cases3+cases4
controls=controls1+controls2+controls3+controls4

# Proportion of cases
s <- cases / (cases + controls)
# Effective sample size
neff <- 4 / (1 / cases + 1 / controls)
# Effective number of cases
neff_cases <- round(neff*s)

# parsing sumstats
#### parsed_input_sumstats <- input_sumstats %>% head()
parsed_input_sumstats <- input_sumstats %>% 
  mutate(
    chromosome=paste0("chr",chromosome),
    position=as.integer(base_pair_location),
    effect_allele=effect_allele,
    non_effect_allele=other_allele,
    effect_size=beta,
    standard_error=standard_error,
    pvalue=p_value,
    sample_size=round(neff),
    n_cases=round(neff_cases)
  ) %>% 
  mutate(initial_variant_id=paste(chromosome,position,non_effect_allele,effect_allele,"b37",sep="_"))
# selecting relevant columns
parsed_input_sumstats <- parsed_input_sumstats %>% select(chromosome,position,initial_variant_id,effect_allele,non_effect_allele,effect_size,standard_error,pvalue,sample_size,n_cases)
parsed_input_sumstats <- na.omit(parsed_input_sumstats)
# save only the first instance of any duplicated initial_variant_id
parsed_input_sumstats <- parsed_input_sumstats[!duplicated(initial_variant_id)]

# preparing liftOver input table
tmp_liftOver <- parsed_input_sumstats %>% mutate(position2=as.integer(position+1)) %>% select(chromosome,position,position2,initial_variant_id)

# writing out and performing liftOver
setwd("/gpfs/data/pierce-lab/james.li/imQTL/tmp/liftOver")
write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
# reading in liftOver mapped results
liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
colnames(liftOver_mapped) <- c("chromosome","position_hg38","position_hg38_2","initial_variant_id")
liftOver_mapped <- liftOver_mapped %>% select(-position_hg38_2)

# joining liftOver coordinates - orientation #1
output_sumstats1 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,non_effect_allele,effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats1 <- output_sumstats1 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# joining liftOver coordinates - orientation #2
output_sumstats2 <- inner_join(liftOver_mapped,(parsed_input_sumstats%>%select(-position)),by=c("chromosome","initial_variant_id")) %>% mutate(panel_variant_id=paste(chromosome,position_hg38,effect_allele,non_effect_allele,"b38",sep="_")) %>% rename(position=position_hg38) %>% mutate(initial_variant_id=panel_variant_id) %>% select(-panel_variant_id) %>% rename(panel_variant_id=initial_variant_id,old_effect=effect_allele,old_non=non_effect_allele) %>% rename(effect_allele=old_non,non_effect_allele=old_effect) %>% mutate(effect_size = -effect_size) %>% mutate(zscore=effect_size/standard_error)
# joining allele frequencies
output_sumstats2 <- output_sumstats2 %>% inner_join(afreq_combined,by=c("panel_variant_id"))

# combining all variants
output_sumstats <- rbind(output_sumstats1,output_sumstats2)

# performing final parsing of sumstats to match the 87 trait GWAS sumstats
output_sumstats <- output_sumstats %>% mutate(
  variant_id=panel_variant_id, 
  current_build="hg38",
  frequency=MAF,
  imputation_status="imputed"
) %>% select(
  variant_id,panel_variant_id,chromosome,position,effect_allele,non_effect_allele,current_build,frequency,sample_size,zscore,pvalue,effect_size,standard_error,imputation_status,n_cases
)

# removing records with na values
output_sumstats <- na.omit(output_sumstats)

# save only the first instance of any duplicated initial_variant_id
output_sumstats <- output_sumstats[!duplicated(panel_variant_id)]

# writing out final output summary statistics
write.table(
  output_sumstats,
  file = gzfile("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/cancer_prostate.txt.gz"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)
