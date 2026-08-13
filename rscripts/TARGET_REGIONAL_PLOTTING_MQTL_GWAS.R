# initializing packages
library(data.table)
library(dplyr)
library(topr)
library(locusplotr)
library(ggplot2)
library(patchwork)
library(ggrepel)

###################################
## PLOTTING CODE: +/- 2000000 BP ##
###################################
set_window_size=500000
half_window_size=set_window_size/2

# specifying imqtl arguments
candidate_str="colon_Epi.cg04093349.chr15:66716105:A:G"
query_rsid="rs56173559"
query_gwas_file="/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/cancer_colon.txt.gz"

# parsing inputs
query_combination=strsplit(candidate_str,split="\\.")[[1]][1]
query_cpg=strsplit(candidate_str,split="\\.")[[1]][2]
query_variant=strsplit(candidate_str,split="\\.")[[1]][3]
query_tissue=strsplit(query_combination, "_")[[1]][1]
parts <- strsplit(query_variant, ":")[[1]]
query_chr <- parts[1]
query_chr_int <- as.integer(gsub("chr","",query_chr))
query_pos <- as.integer(parts[2])

# importing imqtl data
input_filename <- gsub("\\.","_",candidate_str)
input_filename <- paste0(gsub("\\:","_",input_filename),".tsv")
data_imqtl <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_data/",input_filename)) %>% rename(CHROM=CHR,P=pval_gi) %>% select(variant_id,CHROM,POS,P) %>% mutate(Ancestry="imQTL") %>% rename(ID=variant_id)
# importing gwas data
data_gwas <- fread(query_gwas_file)
data_gwas <- data_gwas %>% filter(chromosome==query_chr,position>query_pos-half_window_size,position<query_pos+half_window_size) %>% mutate(variant_id = paste(chromosome,position,non_effect_allele,effect_allele,sep=":"))
data_gwas <- data_gwas %>% rename(CHROM=chromosome,POS=position,P=pvalue) %>% select(variant_id,CHROM,POS,P) %>% mutate(Ancestry="GWAS") %>% rename(ID=variant_id)
data_gwas <- na.omit(data_gwas)
# importing mqtl data
input_filename <- gsub("\\.","_",candidate_str)
input_filename <- paste0(gsub("\\:","_",input_filename),".tsv")
data_mqtl <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_mqtl_data/",input_filename)) %>% rename(CHROM=chr,P=pval_nominal,POS=pos) %>% select(variant_id,CHROM,POS,P) %>% mutate(Ancestry="mQTL") %>% rename(ID=variant_id)

########################################
# COMPUTING LD WITH LEAD VARIANT IMQTL #
########################################
# writing out variant list 
write.table(rbind((data_imqtl%>%select(ID)),data.frame(ID=c(query_variant)))%>%unique(),file="/scratch/jll1/tmp/variant.list",quote=F,row.names=F,col.names=F)
# computing LD
max_snp_window<-nrow(data_imqtl)
system(paste0("module load plink/1.9; plink -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix_bfile --extract /scratch/jll1/tmp/variant.list --ld-snp ",query_variant," --ld-window-kb 500 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
# import pairwise correlations
pairwise_cor <- fread("/scratch/jll1/tmp/tmp_LD.ld")
# Filter the rows where either SNP_A or SNP_B equals query_variant
filtered_df <- pairwise_cor[pairwise_cor$SNP_A == query_variant | pairwise_cor$SNP_B == query_variant, ]
# Create a new column 'OTHER_VARIANT'
filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == query_variant, filtered_df$SNP_B, filtered_df$SNP_A)
# finalizing r2 df
r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(ID=OTHER_VARIANT,r2=R2)
r2_df <- rbind(r2_df,data.frame(ID=query_variant,r2=1))
# joining correlations with main df
tmp_data_imqtl <- left_join(data_imqtl,r2_df,by=c("ID"))
data_imqtl <- tmp_data_imqtl

#######################################
# COMPUTING LD WITH LEAD VARIANT GWAS #
#######################################
# writing out variant list 
write.table(rbind((data_gwas%>%select(ID)),data.frame(ID=c(query_variant)))%>%unique(),file="/scratch/jll1/tmp/variant.list",quote=F,row.names=F,col.names=F)
# computing LD
max_snp_window<-nrow(data_gwas)
system(paste0("module load plink/1.9; plink -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix_bfile --extract /scratch/jll1/tmp/variant.list --ld-snp ",query_variant," --ld-window-kb 500 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
# import pairwise correlations
pairwise_cor <- fread("/scratch/jll1/tmp/tmp_LD.ld")
# Filter the rows where either SNP_A or SNP_B equals query_variant
filtered_df <- pairwise_cor[pairwise_cor$SNP_A == query_variant | pairwise_cor$SNP_B == query_variant, ]
# Create a new column 'OTHER_VARIANT'
filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == query_variant, filtered_df$SNP_B, filtered_df$SNP_A)
# finalizing r2 df
r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(ID=OTHER_VARIANT,r2=R2)
r2_df <- rbind(r2_df,data.frame(ID=query_variant,r2=1))
# joining correlations with main df
tmp_data_gwas <- left_join(data_gwas,r2_df,by=c("ID"))
data_gwas <- tmp_data_gwas

#######################################
# COMPUTING LD WITH LEAD VARIANT MQTL #
#######################################
# writing out variant list 
write.table(rbind((data_mqtl%>%select(ID)),data.frame(ID=c(query_variant)))%>%unique(),file="/scratch/jll1/tmp/variant.list",quote=F,row.names=F,col.names=F)
# computing LD
max_snp_window<-nrow(data_mqtl)
system(paste0("module load plink/1.9; plink -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data/processed_genetic_data_chrprefix_bfile --extract /scratch/jll1/tmp/variant.list --ld-snp ",query_variant," --ld-window-kb 500 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --memory 100000 --out /scratch/jll1/tmp/tmp_LD"))
# import pairwise correlations
pairwise_cor <- fread("/scratch/jll1/tmp/tmp_LD.ld")
# Filter the rows where either SNP_A or SNP_B equals query_variant
filtered_df <- pairwise_cor[pairwise_cor$SNP_A == query_variant | pairwise_cor$SNP_B == query_variant, ]
# Create a new column 'OTHER_VARIANT'
filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == query_variant, filtered_df$SNP_B, filtered_df$SNP_A)
# finalizing r2 df
r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(ID=OTHER_VARIANT,r2=R2)
r2_df <- rbind(r2_df,data.frame(ID=query_variant,r2=1))
# joining correlations with main df
tmp_data_mqtl <- left_join(data_mqtl,r2_df,by=c("ID"))
data_mqtl <- tmp_data_mqtl

######################
# RBIND ALL PLOT DFS #   
######################
plot_df <- rbind(data_imqtl,data_gwas,data_mqtl)

##########################################
# GENERATING COMPREHENSIVE REGIONAL PLOT # 
##########################################
# Ensure Ancestry facet order is fixed
plot_df <- plot_df %>%
  mutate(
    Ancestry = factor(Ancestry, levels = c("GWAS", "mQTL", "imQTL")),
    is_lead = ID == query_variant
  )

# Define x-axis domain (in base pairs)
pos_min <- min(plot_df$POS, na.rm = TRUE)
pos_max <- max(plot_df$POS, na.rm = TRUE)

# Offsets for lead SNP label
x_offset <- 0.02 * ((pos_max - pos_min) / 1e6)
y_offset <- 0.00 * (max(-log10(plot_df$P), na.rm = TRUE))

# Create lead SNP label data
lead_labels_df <- plot_df %>%
  filter(ID == query_variant) %>%
  mutate(rsid = query_rsid) %>%
  group_by(Ancestry) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    x_label = POS / 1e6 + x_offset,
    y_label = -log10(P) + y_offset
  )

current_rsid = query_rsid

# Compute per-panel max -log10(P) and scale 110%
panel_max_df <- plot_df %>%
  group_by(Ancestry) %>%
  summarise(y_max = max(-log10(P), na.rm = TRUE) * 1.10, .groups = "drop") %>%
  mutate(POS = (pos_min + pos_max) / 2,  # middle x
         P = 1,                          # dummy P (so -log10(P)=0)
         dummy = TRUE)

# Add dummy rows to drive y-axis scale
ylim_dummy_df <- panel_max_df %>%
  select(Ancestry, POS, P, y_max) %>%
  mutate(ID = NA, is_lead = FALSE, r2 = NA)

plot_df_final <- bind_rows(plot_df, ylim_dummy_df)

###########################################
library(ggplot2)

p <- ggplot(plot_df_final, aes(x = POS / 1e6, y = -log10(P))) +
  geom_point(data = filter(plot_df_final, !is_lead & !is.na(r2)),
             aes(color = r2), shape = 16, size = 1.8, alpha = 0.8) +
  geom_point(data = filter(plot_df_final, !is_lead & is.na(r2)),
             color = "black", shape = 16, size = 1.8, alpha = 0.8) +
  geom_point(data = filter(plot_df_final, is_lead),
             shape = 18, size = 3, color = "blueviolet") +
  geom_blank(data = panel_max_df, aes(x = POS / 1e6, y = y_max)) +  # forces panel-wise y limit
  facet_wrap(~Ancestry, ncol = 1, scales = "free_y") +
  scale_color_gradientn(
    name = expression(r^2),
    colours = c("darkblue", "lightblue", "chartreuse3", "orange", "red"),
    values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1)),
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 0.27,
      barheight = 3.6,
      direction = "vertical"
    )
  ) +
  scale_y_continuous(
    name = expression(-log[10](italic(P)))
  ) +
  coord_cartesian(
    xlim = c(pos_min, pos_max) / 1e6
  ) +
  labs(x = paste("Position on Chromosome", gsub("chr", "", query_chr), "(in Mb)")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    legend.key.height = unit(0.65, "cm"),
    legend.margin = margin(1, 4, 4, 4),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5)
  ) +
  geom_segment(
    data = lead_labels_df,
    aes(x = x_label, xend = POS / 1e6,
        y = y_label, yend = -log10(P)),
    color = "black", linewidth = 0.3
  ) +
  geom_label(
    data = lead_labels_df,
    aes(x = x_label, y = y_label, label = current_rsid),
    size = 2.5, fontface = "bold",
    fill = "white", color = "black",
    label.size = 0.3, hjust = 0, vjust = 0
  )

# Save the plot
output_path <- paste0(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/target_figure_regional_plots/gwas_mqtl_imqtl/",
  gsub("\\.tsv", "", input_filename),
  ".png"
)

ggsave(output_path, plot = p, width = 3.5, height = 3.75, dpi = 1200)
