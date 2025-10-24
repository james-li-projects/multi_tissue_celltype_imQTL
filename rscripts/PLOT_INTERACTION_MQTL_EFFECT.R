library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

imQTL <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/imQTL/top_assoc/tensorQTL_imQTL_colon_Epi.cis_qtl_top_assoc.txt.gz") %>% filter(pval_adj_bh<0.05)

i=6
current_imQTL = imQTL[i,]
current_imQTL_variant = current_imQTL$variant_id
current_imQTL_cpg = current_imQTL$phenotype_id

########################################
# importing normalized/transformed cell type proportions
celltype <- fread("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/processed_interactions_colon_Epi.txt") %>% mutate(V1=gsub("-",".",V1)) 

########################################
# obtaining genotype counts
system(paste0("echo ",current_imQTL_variant," > /gpfs/data/pierce-lab/james.li/imQTL/tmp/extract_imqtl.list"))

system("plink2 -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/extract_imqtl.list --export Av --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/sample_counts")

genotype_counts <- data.frame(t(fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/sample_counts.traw")))
colnames(genotype_counts)[1] <- "Genotype"
genotype_counts <- genotype_counts %>% mutate(V1=rownames(genotype_counts))

g1 <- paste0(genotype_counts["COUNTED","Genotype"],"/",genotype_counts["COUNTED","Genotype"])
g2 <- paste0(genotype_counts["ALT","Genotype"],"/",genotype_counts["COUNTED","Genotype"])
g3 <- paste0(genotype_counts["ALT","Genotype"],"/",genotype_counts["ALT","Genotype"])

genotype_counts <- genotype_counts[-c(1:6),]
genotype_counts <- genotype_counts %>% mutate(Genotype=as.numeric(Genotype))
genotype_counts <- genotype_counts %>% 
  filter(Genotype %in% c(0,1,2)) %>% 
  mutate(Genotype=ifelse(Genotype==2,g1,Genotype)) %>%
  mutate(Genotype=ifelse(Genotype==1,g2,Genotype)) %>%
  mutate(Genotype=ifelse(Genotype==0,g3,Genotype)) %>% 
  mutate(V1=gsub("0_","",gsub("-",".",V1)))

########################################
# importing methylation data
bed<-fread("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/colon.bed")
t_bed <- data.frame(t(bed %>% filter(phenotype_id==current_imQTL_cpg) %>% select(-"#chr",-start,-end,-phenotype_id)))
colnames(t_bed)[1] <- "mval"
t_bed$V1 <- rownames(t_bed)
t_bed <- t_bed %>% mutate(V1=gsub("-",".",V1))

########################################
# Perform inner joins
combined_df <- t_bed %>%
  inner_join(genotype_counts, by = "V1") %>%
  inner_join(celltype, by = "V1") %>% mutate(`Cell Type Proportion`=ifelse(V2>median(V2),"Upper 50%","Lower 50%"))

# Create the plot
plot <- ggplot(combined_df, aes(x = Genotype, y = mval)) +
  geom_boxplot() +
  facet_wrap(~ `Cell Type Proportion`, ncol = 2) + # Horizontal facets
  theme_classic() +
  labs(
    x = "Genotype",
    y = "DNAm (M-value)",
    title = paste(current_imQTL_variant, current_imQTL_cpg, sep = " - ")
  ) +
  theme(
    axis.text = element_text(size = 14),    # Axis numbers (tick labels)
    axis.title = element_text(size = 16),  # Axis titles
    strip.text = element_text(size = 14),  # Facet labels
    plot.title = element_text(size = 16, face = "bold"), # Plot title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # Add black outline
  )

# Save the plot
output_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/interaction_boxplots/",current_imQTL_variant,"_",current_imQTL_cpg, ".png")
output_path <- gsub("\\:",".",output_path) 
ggsave(output_path, plot, width = 10, height = 6, dpi = 300)

