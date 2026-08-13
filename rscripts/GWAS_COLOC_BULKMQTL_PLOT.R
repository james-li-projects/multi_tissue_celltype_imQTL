library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

######################################
# importing cleaned GWAS trait names #
######################################
info <- fread("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/GWAS_Trait_87_Dictionary/GWASes.87.txt",header=T) %>% mutate(V2_JOIN=V2) %>% mutate(p2=gsub(".txt.gz","",V3)) %>%
  mutate(V2=gsub("Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor: ","",V2)) %>%
  mutate(V2=gsub(" - Only Europeans","",V2)) %>%
  mutate(V2=gsub(" (Discovery and replication cohorts, genders pooled)","",V2)) %>%
  mutate(V2=gsub("Non-cancer illness code, self-reported: ","",V2)) %>%
  mutate(V2=gsub("Vascular/heart problems diagnosed by doctor: ","",V2)) %>%
  mutate(V2=gsub(": Pattern","",V2)) %>%
  mutate(V2=gsub("Diagnoses - main ICD10: G43 ","",V2)) %>%
  mutate(V2=gsub("Diagnoses - main ICD10: G40 ","",V2)) 
info <- info %>%
  mutate(V2=gsub("Vascular/heart problems diagnosed by doctor: ","",V2)) %>%
  mutate(V2=gsub("Attention deficit - Hiperactivity disorder","adhd",V2)) %>%
  mutate(V2=gsub(" \\(chronotype)","",V2)) %>%
  mutate(V2=gsub("Education years \\(Discovery and replication cohorts, genders pooled)","Education years",V2)) %>%
  mutate(V2=gsub("Average number of methylene groups per a double bond","NMR Methylene to Double Bond Ratio",V2))

category_trait <- data.frame(readxl::read_xlsx("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/GWAS_Trait_87_Dictionary/category&trait_cleaned.xlsx")) %>% rename(`Trait Name`=Trait)
cleaned_trait_category_names <- inner_join(info,category_trait,by=c("V2_JOIN"="Trait Name")) %>% select(V2,V3,Category) 


# adding the cancer phenotypes we additionally included
cleaned_trait_category_names <- rbind(
  cleaned_trait_category_names,
  data.frame(V2="Colon cancer",V3="cancer_colon.txt.gz",Category="Cancer"),
  data.frame(V2="Lung cancer",V3="cancer_lung.txt.gz",Category="Cancer"),
  data.frame(V2="Ovarian cancer",V3="cancer_ovary.txt.gz",Category="Cancer"),
  data.frame(V2="Prostate cancer",V3="cancer_prostate.txt.gz",Category="Cancer"),
  data.frame(V2="Leukemia",V3="cancer_leukemia.txt.gz",Category="Cancer")
)

########################################################
########################################################
library(stringr)

# Ensure the column is character
cleaned_trait_category_names$V2 <- as.character(cleaned_trait_category_names$V2)

# Replacement vector
replacements <- c(
  # General replacements
  "counts" = "c\\.",
  "count" = "c\\.",
  "High light scatter" = "HLS",
  "Sum basophil neutrophil" = "Basophil+Neutrophil",
  "Sum eosinophil basophil" = "Eosinophil+Basophil",
  "Sum neutrophil eosinophil" = "Neutrophil+Eosinophil",
  
  # Trait-specific replacements
  "Breast Cancer - Estrogen-receptor positive" = "Breast cancer (ER+)",
  "Breast Cancer - Estrogen-receptor negative" = "Breast cancer (ER-)",
  "Breast Cancer (Overall)" = "Breast cancer (overall)",
  "Coronary Artery Disease" = "Coronary artery disease",
  "Intracraneal Volume" = "Intracraneal volume",
  "Insomnia in both sexes" = "Insomnia",
  "Bone Mineral Density \\(Forearm\\)" = "Bone mineral density",
  "NMR Methylene to Double Bond Ratio" = "CH2:C=C ratio",
  "Birth Weight" = "Birth weight",
  "Total cholesterol in HDL" = "HDL Cholesterol",
  "Triglycerides in IDL" = "IDL Triglycerides",
  "Total cholesterol in LDL" = "LDL Cholesterol",
  "Crohn's Disease" = "Crohn's disease",
  "Inflammatory Bowel Disease" = "Inflammatory bowel",
  "inflammatory bowel disease" = "Inflammatory bowel",
  "Ulcerative Colitis" = "Ulcerative colitis",
  "Standing height" = "Height 1",
  "adhd" = "ADHD",
  "Sleeplessness / insomnia" = "Insomnia",
  "Alzheimer" = "Alzheimer's disease",
  "Systemic Lupus Erythematosus" = "Lupus",
  "Sleep Duration" = "Sleep duration",
  "Fasting Glucose" = "Fasting glucose",
  "Fasting Insulin" = "Fasting insulin",
  "Rheumatoid Arthritis" = "Rheumatoid arthritis",
  "Depressive Symptoms" = "Depressive symptoms",
  "Hayfever, allergic rhinitis or eczema" = "Hayfever+Eczema",
  # Additional replacements (lowercase-insensitive terms should be matched separately if needed)
  "hypertension" = "Hypertension",
  "deep venous thrombosis \\(dvt\\)" = "Deep venous thrombosis",
  "asthma" = "Asthma",
  "irritable bowel syndrome" = "Irritable bowel syndrome",
  "type 1 diabetes" = "Type 1 diabetes",
  "type 2 diabetes" = "Type 2 diabetes",
  "hyperthyroidism/thyrotoxicosis" = "Hyperthyroidism/Thyrotoxicosis",
  "hypothyroidism/myxoedema" = "Hypothyroidism/Myxoedema",
  "psychological/psychiatric problem" = "Psychiatric problems",
  "multiple sclerosis" = "Multiple sclerosis",
  "parkinsons disease" = "Parkinson's disease",
  "migraine" = "Migraine",
  "schizophrenia" = "Schizophrenia",
  "osteoporosis" = "Osteoporosis",
  "ankylosing spondylitis" = "Ankylosing spondylitis",
  "eczema/dermatitis" = "Eczema/Dermatitis",
  "psoriasis" = "Psoriasis",
  "crohns disease" = "Crohn's disease",
  "ulcerative colitis" = "Ulcerative colitis",
  "rheumatoid arthritis" = "Rheumatoid arthritis",
  "gout" = "Gout",
  "high cholesterol" = "High cholesterol",
  "insomnia" = "Insomnia",
  "Body mass index \\(BMI\\)" = "Body mass index",
  "Body fat percentage" = "Body fat perc.",
  "Hair/balding pattern 2" = "Hair loss/Balding",
  "Hair/balding pattern 3" = "Hair loss/Balding",
  "Hair/balding pattern 4" = "Hair loss/Balding",
  "Mother's age at death" = "Mother's age at death",
  "Blood clot in the leg \\(DVT\\)" = "Blood clot (leg)",
  "Blood clot in the lung" = "Blood clot (lung)"
)

# Apply all replacements
cleaned_trait_category_names$V2 <- str_replace_all(cleaned_trait_category_names$V2, replacements)

# special modification for Height 2
cleaned_trait_category_names <- cleaned_trait_category_names %>% mutate(V2 = recode(V2, "Height" = "Height 2"))

# Modify cleaned_trait_category_names directly
cleaned_trait_category_names <- cleaned_trait_category_names %>%
  group_by(V2) %>%
  mutate(dup_count = row_number()) %>%
  ungroup() %>%
  mutate(V2 = ifelse(duplicated(V2) | duplicated(V2, fromLast = TRUE),
                     paste0(V2, " ", dup_count),
                     V2)) %>%
  select(-dup_count)


########################################################

###########################
# importing coloc results #
###########################
# set working directory to coloc output
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/coloc_output")

# obtaining list of coloc output files
coloc_output_list <- list.files()

# importing all coloc output files
gwas_coloc_output_df <- data.frame()
for (coloc_output_file in coloc_output_list) {
  print(coloc_output_file)
  tmp_gwas_coloc_output_df <- fread(coloc_output_file) %>% mutate(V3=gsub("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/","",trait))
  gwas_coloc_output_df <- rbind(gwas_coloc_output_df,tmp_gwas_coloc_output_df)
}

# joining coloc results to trait lists
## joined_gwas_coloc_output_df <- left_join(gwas_coloc_output_df,cleaned_trait_category_names,by=c("V3"))
## gwas_coloc_output_df <- joined_gwas_coloc_output_df
library(data.table)
setDT(gwas_coloc_output_df)
setDT(cleaned_trait_category_names)
setkey(cleaned_trait_category_names, V3)
gwas_coloc_output_df <- cleaned_trait_category_names[
  gwas_coloc_output_df,
  on = "V3"
]

# identifying all mqtls that were included in the coloc analysis
library(data.table)
setDT(gwas_coloc_output_df)
included_mqtl <- unique(
  gwas_coloc_output_df[, .(combination, cpg_id, variant_id)]
)
nrow(included_mqtl)
uniqueN(included_mqtl$variant_id)

# identifying mqtls that exhibited a coloc PP4>0.5 for at least one GWAS trait  
pp4_0.5_mqtl <- gwas_coloc_output_df %>% filter(PP.H4.abf>0.5) %>% rename(Trait=V2) 
unique_pp4_0.5_mqtl <- pp4_0.5_mqtl %>% select(combination, cpg_id, variant_id) %>% unique()
nrow(unique_pp4_0.5_mqtl)
length(unique(unique_pp4_0.5_mqtl$variant_id))

#################################
# creating a unique combination, cpg_id, variant_id list for imQTLs that colocalized with GWAS signals
pp4_0.5_mqtl <- pp4_0.5_mqtl %>% mutate(mqtl_id=paste(combination,cpg_id,variant_id,sep="."))
gwas_coloc_list <- unique((pp4_0.5_mqtl)$mqtl_id)

# identifying mqtls with multiple colocalization signals
counts <- data.frame(pp4_0.5_mqtl %>% select(mqtl_id) %>% table()) %>% arrange(desc(Freq)) 

############################
# plot histogram for number of colocalizing signals
############################

# Set output file path
output_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/GWAS_coloc_multiple_hit_histogram.png"

# Create the histogram plot
png(output_path, width = 5, height = 5, units = "in", res = 300)
ggplot(counts, aes(x = Freq)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(
    x = "Number of colocalizing GWAS traits",
    y = "Count"
  ) +
  theme_classic(base_size = 14)
dev.off()

########################################
library(ggplot2)
library(dplyr)

# Step 1: Count real bars only
pp4_counts <- pp4_0.5_mqtl %>%
  count(Category, Trait) %>%
  mutate(
    Trait = factor(Trait),
    Category = factor(Category)
  )

# Step 2: Plot with facet_grid to simulate grouped spacing
png("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/GWAS_coloc_hist.png", 
    res = 300, unit = "in", height = 4, width = 13)

ggplot(pp4_counts, aes(x = Trait, y = n, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Number of colocalizations") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom", legend.direction = "horizontal",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 30),
    legend.text = element_text(size = 7.5),
    legend.title = element_blank()
  ) + 
  guides(fill = guide_legend(nrow = 1)) + # replace `fill` with your aesthetic
  scale_fill_manual(values = scales::hue_pal()(length(unique(pp4_counts$Category))))
dev.off()


###################################
# simple pie chart of coloc vs not coloc
library(dplyr)
library(ggplot2)

# STEP 1: Aggregate to one row per imQTL (maximum PP.H4.abf per imQTL)
gwas_coloc_output_df <- gwas_coloc_output_df %>% mutate(mqtl_id=paste(combination,cpg_id,variant_id,sep="."))
gwas_coloc_summary <- gwas_coloc_output_df %>%
  group_by(mqtl_id) %>%
  summarise(max_PP.H4 = max(PP.H4.abf, na.rm = TRUE)) %>%
  mutate(status = ifelse(max_PP.H4 > 0.5, "Colocalized with GWAS", "Not colocalized"))

# STEP 2: Prepare data for pie chart
pie_df <- gwas_coloc_summary %>%
  count(status) %>%
  mutate(
    Fraction = n / sum(n),
    Label = paste0(status, "\n", n, " (", round(Fraction * 100, 1), "%)")
  )

# STEP 3: Create pie chart with theme_classic and no title
p <- ggplot(pie_df, aes(x = "", y = n, fill = status)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  ) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("Colocalized with GWAS" = "gold", "Not colocalized" = "gray70"))

# STEP 4: Save to file
ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/coloc_piechart_gwas_mqtl_default_prior.png",
  plot = p,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)


############################################
# Assembling a GWAS coloc table to output 
library(dplyr)
library(tidyr)

# Lookup maps
tissue_map <- c("breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
                "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary")
celltype_map <- c("B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell",
                  "NK" = "NK cell", "Neutro" = "Neutrophil", "Mono" = "Monocyte",
                  "IC" = "Immune cells", "DC" = "Dendritic cell", "Epi" = "Epithelial cell",
                  "Endo" = "Endothelial cell", "EndoC" = "Endothelial cell",
                  "MP" = "Macrophage", "Macro" = "Macrophage", "Fib" = "Fibroblast",
                  "Lym" = "Lymphocyte", "Stromal" = "Stromal cell", "Mye" = "Myeloid cell",
                  "EC" = "Endothelial cell", "SM" = "Smooth muscle cell",
                  "LE" = "Luminal epithelium", "BE" = "Basal epithelium")

#########################################
# WRITING OUT SIGNIFICANT COLOC RESULTS #
#########################################
# Apply the mapping
OUTPUT_TABLE_GWAS_COLOC <- pp4_0.5_mqtl %>% mutate(mqtl_id=paste(combination,cpg_id,variant_id,sep=".")) %>%
  mutate(tissue=combination) %>%
  mutate(
    tissue = tissue_map[tissue]
  ) 
# writing out this table
write.table(OUTPUT_TABLE_GWAS_COLOC,file="/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/OUTPUT_TABLE_GWAS_COLOC.tsv",quote=F,sep="\t",row.names=F,col.names=T)


#################################
# WRITING OUT ALL COLOC RESULTS #
#################################
# Apply the mapping
ALL_OUTPUT_TABLE_GWAS_COLOC <- gwas_coloc_output_df %>% rename(Trait=V2) %>% mutate(mqtl_id=paste(combination,cpg_id,variant_id,sep=".")) %>%
  mutate(tissue=combination) %>%
  mutate(
    tissue = tissue_map[tissue]
  ) 
# writing out this table only if needed [this is gonna be a super large table]
## write.table(ALL_OUTPUT_TABLE_GWAS_COLOC,file="/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/ALL_OUTPUT_TABLE_GWAS_COLOC.tsv",quote=F,sep="\t",row.names=F,col.names=T)










# create histograms of each PP
library(ggplot2)
library(tidyr)
library(dplyr)

# Reshape the data from wide to long for PP.H* columns
pp_long_df <- OUTPUT_TABLE_GWAS_COLOC %>%
  pivot_longer(
    cols = starts_with("PP.H"),
    names_to = "PP_Hypothesis",
    values_to = "Posterior_Probability"
  ) %>%
  mutate(
    PP_Hypothesis = recode(PP_Hypothesis,
                           "PP.H0.abf" = "PP0",
                           "PP.H1.abf" = "PP1",
                           "PP.H2.abf" = "PP2",
                           "PP.H3.abf" = "PP3",
                           "PP.H4.abf" = "PP4")
  )

# Define a color palette for each PP
pp_colors <- c(
  "PP0" = "#1f77b4",  # blue
  "PP1" = "#ff7f0e",  # orange
  "PP2" = "#2ca02c",  # green
  "PP3" = "#d62728",  # red
  "PP4" = "#9467bd"   # purple
)

# Plot faceted histogram with 100 bins and color-coded by PP
pp_hist_plot <- ggplot(pp_long_df, aes(x = Posterior_Probability, fill = PP_Hypothesis)) +
  geom_histogram(bins = 100, color = "black", alpha = 0.7) +
  scale_fill_manual(values = pp_colors) +
  facet_wrap(~ PP_Hypothesis, ncol = 1, scales = "free_y") +
  labs(
    x = "Posterior Probability",
    y = "Count"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

# Save the plot
ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc_bulkmqtl/PP_HIST.png",
  plot = pp_hist_plot,
  width = 6, height = 8, dpi = 300
)
