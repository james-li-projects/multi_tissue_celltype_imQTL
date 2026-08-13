library(data.table)
library(dplyr)





# putative cell type specific mQTLs
candidate_ct_specific_imqtl_DF <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/CANDIDATE_CT_SPECIFIC/tables/candidate_ct_specific_imqtl_DF.tsv")
putativeCTS_imQTL <- unique(candidate_ct_specific_imqtl_DF$target_cpg)

# joining cell type specific coloc results with directionality and parsed effects
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")
thirdclass_imQTL <- unique((wide_parsed_imqtl%>%filter(Category=="Unknown"))$phenotype_id)

# only imQTL, not bulk mQTL
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/imQTL_cpg_df.RData")
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_cpg_df.RData")
only_imQTL <- setdiff(unique(imQTL_cpg_df$phenotype_id),unique(mQTL_cpg_df$phenotype_id))

length(putativeCTS_imQTL)
length(thirdclass_imQTL)
length(only_imQTL)


# install.packages("eulerr")
library(eulerr)

# Make sure vectors are unique
putativeCTS_imQTL <- unique(putativeCTS_imQTL)
thirdclass_imQTL  <- unique(thirdclass_imQTL)
only_imQTL        <- unique(only_imQTL)

# Create Euler/Venn object
fit <- euler(list(
  "Putative CT-specific\nmQTLs" = putativeCTS_imQTL,
  "Identified as imQTL, no detectable\nmarginal mQTL" = thirdclass_imQTL
))

# High-resolution PNG
png(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/CONSOLIDATE_CELLTYPESPEC/venn.png",
  width = 3800,
  height = 3000,
  res = 400
)

plot(
  fit,
  quantities = list(
    type = "counts",
    fontsize = 15
  ),
  labels = list(
    fontsize = 15
  ),
  fills = list(
    fill = c("#4E79A7", "#E15759"),
    alpha = 0.55
  ),
  edges = list(
    col = "black",
    lwd = 2
  ),
  main = ""
)

dev.off()

