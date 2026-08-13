library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

###################################
# defining input and output paths #
###################################
# input query mqtl chunk containing full sumstats for a chromosome
## example: current_qtl_chunk_filename <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/mQTL/full_assoc_tsv/tensorQTL_mQTL_prostate.cis_qtl_pairs.chr22.tsv"
args <- commandArgs(trailingOnly = TRUE)
current_qtl_chunk_filename <- args[1]
current_qtl_chunk <- fread(current_qtl_chunk_filename)
# output mqtl window directory
out_dir_mqtl <- "/scratch/jll1/imQTL/coloc_windows/extracted_mqtl_windows"
dir.create(out_dir_mqtl, recursive = TRUE, showWarnings = FALSE)
# output gwas window directory
out_dir_gwas <- "/scratch/jll1/imQTL/coloc_windows/extracted_gwas_windows"
dir.create(out_dir_gwas, recursive = TRUE, showWarnings = FALSE)

# extracting current dataset, chromosome, tissue, and tissue-specific top association filename
current_dataset <- basename(dirname(dirname(dirname(current_qtl_chunk_filename))))
current_chr <- sub(".*\\.chr([0-9XYM]+)\\.tsv$", "\\1", basename(current_qtl_chunk_filename))
current_tissue <- sub(".*tensorQTL_mQTL_([^.]*)\\..*", "\\1", current_qtl_chunk_filename)
qtl_top_assoc_filename <- file.path(
  dirname(dirname(current_qtl_chunk_filename)),
  "top_assoc",
  sub("\\.cis_qtl_pairs\\.chr[0-9XYM]+\\.tsv$", ".cis_qtl.txt.gz",
      basename(current_qtl_chunk_filename))
)
paste(current_dataset,current_chr,qtl_top_assoc_filename)

# importing the top association file 
current_qtl_top_assoc <- fread(qtl_top_assoc_filename) %>% mutate(pval_adj_bh=p.adjust(pval_beta,method="fdr")) %>% filter(pval_adj_bh < 0.05)

########################
# Extract mQTL windows #
########################
# parsing tables for extracting windows
library(data.table)
setDT(current_qtl_chunk)
setDT(current_qtl_top_assoc)
current_qtl_chunk[
  , c("variant_chr", "variant_pos") :=
    tstrsplit(variant_id, ":", fixed = TRUE, keep = 1:2)
][
  , variant_pos := as.integer(variant_pos)
]

current_qtl_top_assoc[
  , c("variant_chr", "lead_pos") :=
    tstrsplit(variant_id, ":", fixed = TRUE, keep = 1:2)
][
  , lead_pos := as.integer(lead_pos)
]

# Only retain top associations from the chromosome represented by this chunk
current_qtl_top_assoc_chr <- current_qtl_top_assoc[
  variant_chr == paste0("chr", current_chr)
]

setkey(current_qtl_chunk, phenotype_id, variant_chr)

for (i in seq_len(nrow(current_qtl_top_assoc_chr))) {
  
  phenotype    <- current_qtl_top_assoc_chr$phenotype_id[i]
  chromosome   <- current_qtl_top_assoc_chr$variant_chr[i]
  lead_pos     <- current_qtl_top_assoc_chr$lead_pos[i]
  lead_variant <- current_qtl_top_assoc_chr$variant_id[i]
  
  window_dt <- current_qtl_chunk[
    .(phenotype, chromosome),
    on = .(phenotype_id, variant_chr),
    nomatch = 0L
  ][
    variant_pos >= lead_pos - 250000L &
      variant_pos <= lead_pos + 250000L
  ]
  
  fwrite(
    window_dt[, !c("variant_chr", "variant_pos")],
    file.path(
      out_dir_mqtl,
      sprintf(
        "%s_%s_%s.tsv",
        current_tissue,
        phenotype,
        gsub(":", "_", lead_variant, fixed = TRUE)
      )
    ),
    sep = "\t"
  )
}
