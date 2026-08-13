library(data.table)
library(dplyr)
imQTL <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/OUTPUT_TABLE_GWAS_COLOC.tsv")
mQTL <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/OUTPUT_TABLE_GWAS_COLOC.tsv")

nrow(imQTL)
unique_imQTL_cpg <- unique(imQTL$cpg_id)
length(unique_imQTL_cpg)

nrow(mQTL)
unique_mQTL_cpg <- unique(mQTL$cpg_id)
length(unique_mQTL_cpg)

length(intersect(unique_imQTL_cpg,unique_mQTL_cpg))
length(setdiff(unique_imQTL_cpg,unique_mQTL_cpg))
length(setdiff(unique_mQTL_cpg,unique_imQTL_cpg))

# All categories appearing in either dataset
categories <- union(
  unique(as.character(imQTL$Category)),
  unique(as.character(mQTL$Category))
)

# Remove missing categories, if present
categories <- categories[!is.na(categories)]

# Total numbers of imQTL and mQTL observations
n_imQTL <- nrow(imQTL)
n_mQTL  <- nrow(mQTL)

# Run Fisher's exact test for each category
fisher_results <- do.call(
  rbind,
  lapply(categories, function(category_i) {
    
    imQTL_in  <- sum(imQTL$Category == category_i, na.rm = TRUE)
    imQTL_out <- n_imQTL - imQTL_in
    
    mQTL_in   <- sum(mQTL$Category == category_i, na.rm = TRUE)
    mQTL_out  <- n_mQTL - mQTL_in
    
    contingency_table <- matrix(
      c(
        imQTL_in,  imQTL_out,
        mQTL_in,   mQTL_out
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        QTL_type = c("imQTL", "mQTL"),
        Category_status = c("In_category", "Not_in_category")
      )
    )
    
    fisher_test <- fisher.test(contingency_table)
    
    data.frame(
      Category = category_i,
      imQTL_in_category = imQTL_in,
      imQTL_not_in_category = imQTL_out,
      mQTL_in_category = mQTL_in,
      mQTL_not_in_category = mQTL_out,
      imQTL_proportion = imQTL_in / n_imQTL,
      mQTL_proportion = mQTL_in / n_mQTL,
      odds_ratio = unname(fisher_test$estimate),
      conf_low = fisher_test$conf.int[1],
      conf_high = fisher_test$conf.int[2],
      p_value = fisher_test$p.value,
      stringsAsFactors = FALSE
    )
  })
)

# Multiple-testing correction
fisher_results$FDR <- p.adjust(
  fisher_results$p_value,
  method = "BH"
)

# Direction of enrichment
fisher_results$enrichment <- ifelse(
  fisher_results$odds_ratio > 1,
  "Enriched in imQTL",
  ifelse(
    fisher_results$odds_ratio < 1,
    "Depleted in imQTL",
    "No difference"
  )
)

# Sort by FDR
fisher_results <- fisher_results[
  order(fisher_results$FDR),
]

rownames(fisher_results) <- NULL

fisher_results %>% arrange(desc(imQTL_proportion)) %>% select(Category, imQTL_proportion, mQTL_proportion, p_value, FDR, enrichment) %>% filter(enrichment=="Enriched in imQTL")

sort(table(imQTL$Category))
sort(table(mQTL$Category))

write.table(fisher_results,file="/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/compare_imqtl_mqtl_gwascategory.tsv",quote=F,sep="\t",row.names=F,col.names=T)
