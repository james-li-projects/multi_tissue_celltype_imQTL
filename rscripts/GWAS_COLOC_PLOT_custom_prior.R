library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

##################################
# importing eqtl colocalizations #
##################################
all_eqtl_imqtl_coloc <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/all_eqtl_imqtl_coloc.tsv") %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
eqtl_coloc_list <- (all_eqtl_imqtl_coloc %>% filter(coloc_0.5_all==TRUE))$imqtl_id

######################################
# importing cleaned GWAS trait names #
######################################
info <- fread("/gpfs/data/huo-lab/BCAC/james.li/GWAS_Trait_87_Dictionary/GWASes.87.txt",header=F) %>% mutate(V2_JOIN=V2) %>% mutate(p2=gsub(".txt.gz","",V3)) %>%
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

category_trait <- data.frame(readxl::read_xlsx("/gpfs/data/huo-lab/BCAC/james.li/GWAS_Trait_87_Dictionary/category&trait_cleaned.xlsx")) %>% rename(`Trait Name`=Trait)
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
  "Hayfever, allergic rhinitis or eczema" = "Allergies (Hayfever/Eczema)",
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
########################################################




###########################
# loading cpgs mapped to mqtls
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_cpg_df.RData")
mQTL_cpg_df <- mQTL_cpg_df %>%
  mutate(
    Tissue = case_when(
      Tissue == "Colon" ~ "colon",
      Tissue == "Lung" ~ "lung",
      Tissue == "Ovary" ~ "ovary",
      Tissue == "Prostate" ~ "prostate",
      Tissue == "Whole Blood" ~ "wb",
      TRUE ~ Tissue
    )
  ) %>% mutate(tissue_cpg=paste(Tissue,phenotype_id,sep="_"))

###########################
# importing coloc results #
###########################
# set working directory to coloc output
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/coloc_output_custom_prior")

# obtaining list of coloc output files
coloc_output_list <- list.files()

# importing all coloc output files
gwas_coloc_output_df <- data.frame()
for (coloc_output_file in coloc_output_list) {
  print(coloc_output_file)
  tmp_gwas_coloc_output_df <- fread(coloc_output_file) %>% mutate(V3=gsub("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/","",trait))
  gwas_coloc_output_df <- rbind(gwas_coloc_output_df,tmp_gwas_coloc_output_df)
}

# identifying which coloc tissue/cpgs were mapped as mQTLs
gwas_coloc_output_df <- gwas_coloc_output_df %>% separate(combination,into=c("tissue","celltype"),remove=F,sep="_") %>% mutate(tissue_cpg=paste(tissue,cpg_id,sep="_"))
gwas_coloc_output_df <- gwas_coloc_output_df %>% mutate(
  mqtl_status = ifelse(tissue_cpg %in% mQTL_cpg_df$tissue_cpg,"Yes","No")
)

# joining coloc results to trait lists
joined_gwas_coloc_output_df <- left_join(gwas_coloc_output_df,cleaned_trait_category_names,by=c("V3"))
gwas_coloc_output_df <- joined_gwas_coloc_output_df

# identifying all imqtls that were included in the coloc analysis
included_imqtl <- gwas_coloc_output_df %>% select(combination, cpg_id, variant_id) %>% unique()
length(unique(included_imqtl$variant_id))

# identifying imqtls that exhibited a coloc PP4>0.5 for at least one GWAS trait  
pp4_0.5_imqtl <- gwas_coloc_output_df %>% filter(PP.H4.abf>0.5) %>% rename(Trait=V2) 
unique_pp4_0.5_imqtl <- pp4_0.5_imqtl %>% select(combination, cpg_id, variant_id) %>% unique()
length(unique(unique_pp4_0.5_imqtl$variant_id))

#################################
# creating a unique combination, cpg_id, variant_id list for imQTLs that colocalized with GWAS signals
pp4_0.5_imqtl <- pp4_0.5_imqtl %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
gwas_coloc_list <- (pp4_0.5_imqtl)$imqtl_id

# identifying imqtls with multiple colocalization signals
counts <- data.frame(pp4_0.5_imqtl %>% select(imqtl_id) %>% table()) %>% arrange(desc(Freq)) 

############################
# plot violin plot for number of colocalizing signals
############################

# Add a dummy x variable for violin grouping
counts$group <- ""

# Set output file path
output_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/GWAS_coloc_violin.png"

# Create the violin plot
png(output_path, width = 5, height = 5, units = "in", res = 300)
ggplot(counts, aes(x = group, y = Freq)) +
  geom_violin(fill = "skyblue", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.5) +
  labs(
    y = "Number of colocalizing GWAS traits",
    x = ""
  ) +
  theme_classic(base_size = 14) + 
  theme(axis.ticks.x = element_blank())
dev.off()

########################################
library(ggplot2)
library(dplyr)

# Step 1: Count real bars only
pp4_counts <- pp4_0.5_imqtl %>%
  count(Category, Trait) %>%
  mutate(
    Trait = factor(Trait),
    Category = factor(Category)
  )

# Step 2: Plot with facet_grid to simulate grouped spacing
png("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/GWAS_coloc_hist_custom_prior.png", 
    res = 300, unit = "in", height = 5, width = 20)

ggplot(pp4_counts, aes(x = Trait, y = n, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Number of colocalizations") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines"),
    legend.position = "bottom", legend.direction = "horizontal",
    plot.margin = margin(t = 10, r = 30, b = 10, l = 30)
  ) + 
  guides(fill = guide_legend(nrow = 1)) + # replace `fill` with your aesthetic
  scale_fill_manual(values = scales::hue_pal()(length(unique(pp4_counts$Category))))

dev.off()


##############################
# Load required libraries
library(dplyr)
library(ggplot2)

# STEP 1: Aggregate to one row per imQTL (maximum PP.H4.abf per imQTL)
gwas_coloc_output_df <- gwas_coloc_output_df %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
gwas_coloc_summary <- gwas_coloc_output_df %>%
  group_by(imqtl_id) %>%
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
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/coloc_piechart_gwas_imqtl_custom_prior.png",
  plot = p,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)
