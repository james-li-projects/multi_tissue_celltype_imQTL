library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#######################################
# paths
#######################################

imqtl_top_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc"

jobs_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/SCEQTL/40154479_Wang"

blueprint_neutro_file <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/SCEQTL/blueprint/QTD000026.all.tsv.gz"

out_file <- "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/imqtl_variant_sceqtl_summary_stats_wide.tsv"

#######################################
# 1. import significant imQTL lead variants
#######################################

file_list <- list.files(
  imqtl_top_dir,
  pattern = "\\.cis_qtl_top_assoc\\.txt\\.gz$",
  full.names = TRUE
)

imqtl_list <- list()

for (f in file_list) {
  
  message("Reading imQTL file: ", basename(f))
  
  tmp <- fread(f)
  
  tmp <- tmp %>%
    mutate(
      pval_adj_bonf = p.adjust(pval_emt, method = "bonferroni"),
      combination = basename(f),
      Dataset = "HEALS"
    ) %>%
    filter(pval_adj_bonf < 0.05) %>%
    select(
      imqtl_cpg_id = phenotype_id,
      imqtl_variant_id = variant_id,
      imqtl_pval_emt = pval_emt,
      imqtl_pval_adj_bonf = pval_adj_bonf,
      combination,
      Dataset
    )
  
  imqtl_list[[basename(f)]] <- tmp
}

imqtl_leads <- rbindlist(imqtl_list, fill = TRUE)

#######################################
# 2. parse imQTL variant IDs and cell type
#######################################

imqtl_leads <- imqtl_leads %>%
  separate(
    imqtl_variant_id,
    into = c("chr", "pos", "imqtl_a2", "imqtl_a1"),
    sep = ":",
    remove = FALSE
  ) %>%
  mutate(
    pos = as.integer(pos),
    chunk_num = as.integer(gsub("chr", "", chr)),
    imqtl_ct = gsub(
      "tensorQTL_imQTL_", "",
      gsub("\\.cis_qtl_top_assoc\\.txt\\.gz$", "", combination)
    )
  ) %>%
  filter(imqtl_ct != "wb_Eosino") %>%
  distinct()

setDT(imqtl_leads)

#######################################
# 3. make imQTL variant lookup
#######################################
# This allows matching either allele order:
# chr:pos:a2:a1 or chr:pos:a1:a2

imqtl_lookup <- imqtl_leads[
  ,
  .(
    imqtl_cpg_id,
    imqtl_variant_id,
    imqtl_ct,
    imqtl_chr = chr,
    imqtl_pos = pos,
    imqtl_a1,
    imqtl_a2,
    imqtl_pval_emt,
    imqtl_pval_adj_bonf,
    match_id_forward = paste(chr, pos, imqtl_a2, imqtl_a1, sep = ":"),
    match_id_reverse = paste(chr, pos, imqtl_a1, imqtl_a2, sep = ":")
  )
]

lookup_forward <- imqtl_lookup[
  ,
  .(
    imqtl_cpg_id,
    imqtl_variant_id,
    imqtl_ct,
    imqtl_chr,
    imqtl_pos,
    imqtl_a1,
    imqtl_a2,
    imqtl_pval_emt,
    imqtl_pval_adj_bonf,
    match_variant_id = match_id_forward,
    allele_match = "forward"
  )
]

lookup_reverse <- imqtl_lookup[
  ,
  .(
    imqtl_cpg_id,
    imqtl_variant_id,
    imqtl_ct,
    imqtl_chr,
    imqtl_pos,
    imqtl_a1,
    imqtl_a2,
    imqtl_pval_emt,
    imqtl_pval_adj_bonf,
    match_variant_id = match_id_reverse,
    allele_match = "reverse"
  )
]

variant_lookup <- rbindlist(list(lookup_forward, lookup_reverse), fill = TRUE)
setkey(variant_lookup, match_variant_id)

#######################################
# 4. import Wang/JOBS sceQTL files
#######################################

jobs_files <- list.files(
  jobs_dir,
  pattern = "_JOBS_summary_statistics\\.txt\\.gz$",
  full.names = TRUE
)

eqtl_match_list <- list()

for (f in jobs_files) {
  
  current_sceqtl_ct <- gsub("_JOBS_summary_statistics\\.txt\\.gz$", "", basename(f))
  message("Reading JOBS sceQTL file: ", current_sceqtl_ct)
  
  eqtl <- fread(f)
  
  eqtl_tmp <- eqtl[
    ,
    .(
      eqtl_gene_id = gene,
      eqtl_variant_id = paste0(
        "chr",
        gsub("_", ":", sub("_b38$", "", snp_hg38))
      ),
      sceqtl_ct = current_sceqtl_ct,
      eqtl_beta = beta.est,
      eqtl_se = beta.se,
      eqtl_pval = p_val
    )
  ]
  
  eqtl_tmp[, eqtl_varbeta := eqtl_se^2]
  
  setDT(eqtl_tmp)
  setkey(eqtl_tmp, eqtl_variant_id)
  
  matched <- merge(
    variant_lookup,
    eqtl_tmp,
    by.x = "match_variant_id",
    by.y = "eqtl_variant_id",
    all = FALSE,
    allow.cartesian = TRUE
  )
  
  eqtl_match_list[[current_sceqtl_ct]] <- matched
}

#######################################
# 5. import Blueprint neutrophil sceQTL file
#######################################

message("Reading Blueprint neutrophil sceQTL file")

blueprint_eqtl <- fread(blueprint_neutro_file)

blueprint_tmp <- blueprint_eqtl[
  ,
  .(
    eqtl_gene_id = molecular_trait_id,
    eqtl_variant_id = chartr("_", ":", variant),
    sceqtl_ct = "Neutrophil_Blueprint",
    eqtl_beta = beta,
    eqtl_se = se,
    eqtl_pval = pvalue
  )
]

blueprint_tmp[, eqtl_varbeta := eqtl_se^2]

setDT(blueprint_tmp)
setkey(blueprint_tmp, eqtl_variant_id)

blueprint_matched <- merge(
  variant_lookup,
  blueprint_tmp,
  by.x = "match_variant_id",
  by.y = "eqtl_variant_id",
  all = FALSE,
  allow.cartesian = TRUE
)

eqtl_match_list[["Neutrophil_Blueprint"]] <- blueprint_matched

#######################################
# 6. combine long-format matched eQTLs
#######################################

eqtl_long <- rbindlist(eqtl_match_list, fill = TRUE)

#######################################
# 7. harmonize beta direction
#######################################
# Raw beta is kept.
# Harmonized beta flips sign when the eQTL allele order is reversed.

eqtl_long[
  ,
  eqtl_beta_harmonized := fifelse(
    allele_match == "reverse",
    -eqtl_beta,
    eqtl_beta
  )
]

#######################################
# 8. pivot to massively wide dataframe
#######################################

eqtl_wide <- dcast(
  eqtl_long,
  imqtl_cpg_id +
    imqtl_variant_id +
    imqtl_ct +
    imqtl_chr +
    imqtl_pos +
    imqtl_a1 +
    imqtl_a2 +
    imqtl_pval_emt +
    imqtl_pval_adj_bonf +
    eqtl_gene_id ~ sceqtl_ct,
  value.var = c(
    "eqtl_beta",
    "eqtl_beta_harmonized",
    "eqtl_se",
    "eqtl_varbeta",
    "eqtl_pval"
  )
)

#######################################
# 9. save output
#######################################

fwrite(eqtl_wide, out_file, sep = "\t")

message("Saved wide dataframe to: ", out_file)
message("Rows: ", nrow(eqtl_wide))
message("Columns: ", ncol(eqtl_wide))


# removing rows with NA values
eqtl_wide_noNA <- na.omit(eqtl_wide)
head(eqtl_wide_noNA,1)


#######################################
# 10. identifying cell type specific records
#######################################


library(data.table)
library(dplyr)

setDT(eqtl_wide_noNA)

#######################################
# thresholds
#######################################

nominal_p_cutoff <- 0.05
specific_p_cutoff <- 1e-3

#######################################
# p-value columns
#######################################

pval_cols <- grep("^eqtl_pval_", names(eqtl_wide_noNA), value = TRUE)

cell_types <- gsub("^eqtl_pval_", "", pval_cols)

#######################################
# define broad cell type groups
#######################################

cell_type_group_map <- data.table(
  cell_type = cell_types,
  broad_cell_type = case_when(
    grepl("^B_", cell_types) ~ "B_cell",
    grepl("^CD4_", cell_types) ~ "CD4_T_cell",
    grepl("^CD8_", cell_types) ~ "CD8_T_cell",
    grepl("^Mono_", cell_types) ~ "Monocyte",
    grepl("^NK", cell_types) ~ "NK_cell",
    grepl("^DC", cell_types) ~ "Dendritic_cell",
    grepl("^Plasma", cell_types) ~ "Plasma_cell",
    grepl("^Neutrophil", cell_types) ~ "Neutrophil",
    TRUE ~ cell_types
  )
)

#######################################
# reshape p-values to long format
#######################################

pval_long <- melt(
  eqtl_wide_noNA,
  id.vars = c(
    "imqtl_cpg_id",
    "imqtl_variant_id",
    "imqtl_ct",
    "imqtl_chr",
    "imqtl_pos",
    "imqtl_a1",
    "imqtl_a2",
    "imqtl_pval_emt",
    "imqtl_pval_adj_bonf",
    "eqtl_gene_id"
  ),
  measure.vars = pval_cols,
  variable.name = "pval_col",
  value.name = "eqtl_pval"
)

pval_long[
  ,
  cell_type := gsub("^eqtl_pval_", "", pval_col)
]

pval_long <- merge(
  pval_long,
  cell_type_group_map,
  by = "cell_type",
  all.x = TRUE
)

#######################################
# identify specificity by broad group
#######################################

specificity_summary <- pval_long[
  ,
  .(
    min_pval = min(eqtl_pval, na.rm = TRUE),
    best_cell_type = cell_type[which.min(eqtl_pval)],
    best_broad_cell_type = broad_cell_type[which.min(eqtl_pval)],
    n_nominal_sig_cell_types = sum(eqtl_pval < nominal_p_cutoff),
    nominal_sig_cell_types = paste(cell_type[eqtl_pval < nominal_p_cutoff], collapse = ";"),
    nominal_sig_broad_cell_types = paste(unique(broad_cell_type[eqtl_pval < nominal_p_cutoff]), collapse = ";"),
    n_nominal_sig_broad_cell_types = uniqueN(broad_cell_type[eqtl_pval < nominal_p_cutoff])
  ),
  by = .(
    imqtl_cpg_id,
    imqtl_variant_id,
    imqtl_ct,
    imqtl_chr,
    imqtl_pos,
    imqtl_a1,
    imqtl_a2,
    imqtl_pval_emt,
    imqtl_pval_adj_bonf,
    eqtl_gene_id
  )
]

specificity_summary[
  ,
  cell_type_specific := (
    min_pval < specific_p_cutoff &
      n_nominal_sig_broad_cell_types == 1
  )
]

#######################################
# keep only cell-type-specific rows
#######################################

cell_type_specific_rows <- specificity_summary[
  cell_type_specific == TRUE
]

#######################################
# optional: merge annotation back onto wide dataframe
#######################################

eqtl_wide_specific_annotated <- merge(
  eqtl_wide_noNA,
  specificity_summary[
    ,
    .(
      imqtl_cpg_id,
      imqtl_variant_id,
      imqtl_ct,
      eqtl_gene_id,
      cell_type_specific,
      best_cell_type,
      best_broad_cell_type,
      min_pval,
      nominal_sig_cell_types,
      nominal_sig_broad_cell_types
    )
  ],
  by = c(
    "imqtl_cpg_id",
    "imqtl_variant_id",
    "imqtl_ct",
    "eqtl_gene_id"
  ),
  all.x = TRUE
)

eqtl_wide_specific_only <- eqtl_wide_specific_annotated[
  cell_type_specific == TRUE
]

table(eqtl_wide_specific_only$best_cell_type)
write.table(eqtl_wide_specific_only, file="/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/ct_spec_eqtl.tsv",quote=F,col.names=T,row.names=F,sep="\t")


#############################################################
#############################################################
# Importing coloc data for all sceQTLs and imQTLs (default priors) #
#############################################################
#############################################################
library(data.table)
library(dplyr)
library(tidyr)

# reimporting ct specific eQTLs
eqtl_wide_specific_only <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/ct_spec_eqtl.tsv")

# Define the directory path
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/default_prior"

# List coloc output files only
file_list <- list.files(
  input_dir,
  pattern = "^coloc_output_.*\\.tsv$",
  full.names = TRUE
)

# Read all files with fread and combine into one data.table
coloc_data <- rbindlist(lapply(file_list, fread), fill = TRUE)

# Create coloc_0.5_all indicator for each sc-eQTL coloc row
coloc_data[, coloc_0.5_all := PP.H4.abf > 0.5]

# creating pair_id variable
coloc_data <- coloc_data %>%
  mutate(pair_id = paste(cpg_id, imqtl_variant_id, sep = "_"))

# separate parsed_combination from tissue and celltype
coloc_data <- coloc_data %>%
  separate(parsed_combination, into = c("tissue", "celltype"), sep = "_", remove = FALSE)

# View result
print(head(coloc_data))


#############################################################
# Filter coloc_data to cell-type-specific sc-eQTLs only
#############################################################

specific_eqtl_filter <- eqtl_wide_specific_only %>%
  filter(cell_type_specific == TRUE) %>%
  transmute(
    cpg_id = imqtl_cpg_id,
    imqtl_variant_id = imqtl_variant_id,
    parsed_combination = imqtl_ct,
    eqtl_phenotype_id = eqtl_gene_id,
    sceqtl_ct = best_cell_type
  ) %>%
  distinct()

# Keep only coloc rows matching cell-type-specific eQTLs
ct_spec_coloc_data <- coloc_data %>%
  inner_join(
    specific_eqtl_filter,
    by = c(
      "cpg_id",
      "imqtl_variant_id",
      "parsed_combination",
      "eqtl_phenotype_id",
      "sceqtl_ct"
    )
  )

print(dim(ct_spec_coloc_data))

ct_spec_coloc_data_pp40.5<-ct_spec_coloc_data %>% filter(PP.H4.abf>0.5)%>%rename(combination=parsed_combination,phenotype_id=cpg_id,variant_id=imqtl_variant_id)

# joining cell type specific coloc results with directionality and parsed effects
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")

joined_ct_spec_coloc_df <- ct_spec_coloc_data_pp40.5 %>% inner_join(wide_parsed_imqtl,by=c("variant_id","phenotype_id","combination","tissue","celltype"))

# saving a tsv of coloc results for cell type specific eQTLs
write.table(joined_ct_spec_coloc_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/default_prior/CT_SPECIFIC_COLOC_OUTPUT_default_prior.tsv",quote=F,row.names=F,col.names=T,sep="\t")



#############################################################
#############################################################
# Importing coloc data for all sceQTLs and imQTLs (custom priors) #
#############################################################
#############################################################
library(data.table)
library(dplyr)
library(tidyr)

# reimporting ct specific eQTLs
eqtl_wide_specific_only <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/ct_spec_eqtl.tsv")

# Define the directory path
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/custom_prior"

# List coloc output files only
file_list <- list.files(
  input_dir,
  pattern = "^coloc_output_.*\\.tsv$",
  full.names = TRUE
)

# Read all files with fread and combine into one data.table
coloc_data <- rbindlist(lapply(file_list, fread), fill = TRUE)

# Create coloc_0.5_all indicator for each sc-eQTL coloc row
coloc_data[, coloc_0.5_all := PP.H4.abf > 0.5]

# creating pair_id variable
coloc_data <- coloc_data %>%
  mutate(pair_id = paste(cpg_id, imqtl_variant_id, sep = "_"))

# separate parsed_combination from tissue and celltype
coloc_data <- coloc_data %>%
  separate(parsed_combination, into = c("tissue", "celltype"), sep = "_", remove = FALSE)

# View result
print(head(coloc_data))


#############################################################
# Filter coloc_data to cell-type-specific sc-eQTLs only
#############################################################

specific_eqtl_filter <- eqtl_wide_specific_only %>%
  filter(cell_type_specific == TRUE) %>%
  transmute(
    cpg_id = imqtl_cpg_id,
    imqtl_variant_id = imqtl_variant_id,
    parsed_combination = imqtl_ct,
    eqtl_phenotype_id = eqtl_gene_id,
    sceqtl_ct = best_cell_type
  ) %>%
  distinct()

# Keep only coloc rows matching cell-type-specific eQTLs
ct_spec_coloc_data <- coloc_data %>%
  inner_join(
    specific_eqtl_filter,
    by = c(
      "cpg_id",
      "imqtl_variant_id",
      "parsed_combination",
      "eqtl_phenotype_id",
      "sceqtl_ct"
    )
  )

print(dim(ct_spec_coloc_data))

ct_spec_coloc_data_pp40.5<-ct_spec_coloc_data %>% filter(PP.H4.abf>0.5)%>%rename(combination=parsed_combination,phenotype_id=cpg_id,variant_id=imqtl_variant_id)

# joining cell type specific coloc results with directionality and parsed effects
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")

joined_ct_spec_coloc_df <- ct_spec_coloc_data_pp40.5 %>% inner_join(wide_parsed_imqtl,by=c("variant_id","phenotype_id","combination","tissue","celltype"))

# saving a tsv of coloc results for cell type specific eQTLs
write.table(joined_ct_spec_coloc_df,file="/gpfs/data/pierce-lab/james.li/imQTL/output/SCEQTL_coloc/custom_prior/CT_SPECIFIC_COLOC_OUTPUT_custom_prior.tsv",quote=F,row.names=F,col.names=T,sep="\t")

