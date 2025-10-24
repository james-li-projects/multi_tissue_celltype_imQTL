library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

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


#################################
# importing imqtl coloc results #
#################################
# set working directory to coloc output
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/coloc_output")

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
joined_gwas_coloc_output_df <- left_join(gwas_coloc_output_df,cleaned_trait_category_names,by=c("V3"))
gwas_coloc_output_df <- joined_gwas_coloc_output_df

# identifying all imqtls that were included in the coloc analysis
included_imqtl <- gwas_coloc_output_df %>% select(combination, cpg_id, variant_id) %>% unique()
nrow(included_imqtl)
length(unique(included_imqtl$variant_id))

# identifying imqtls that exhibited a coloc PP4>0.5 for at least one GWAS trait  
pp4_0.5_imqtl <- gwas_coloc_output_df %>% filter(PP.H4.abf>0.5) %>% rename(Trait=V2) %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
unique_pp4_0.5_imqtl <- pp4_0.5_imqtl %>% select(combination, cpg_id, variant_id) %>% unique()
nrow(unique_pp4_0.5_imqtl)
length(unique(unique_pp4_0.5_imqtl$variant_id))



################################
# importing mqtl coloc results #
################################
# set working directory to coloc output
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/coloc_output_mqtl")

# obtaining list of coloc output files
coloc_output_mqtl_list <- list.files()

# importing all coloc output files
gwas_coloc_output_mqtl_df <- data.frame()
for (coloc_output_mqtl_file in coloc_output_mqtl_list) {
  print(coloc_output_mqtl_file)
  tmp_gwas_coloc_output_mqtl_df <- fread(coloc_output_mqtl_file) %>% mutate(V3=gsub("/gpfs/data/pierce-lab/james.li/GWAS_Trait_87/","",trait))
  gwas_coloc_output_mqtl_df <- rbind(gwas_coloc_output_mqtl_df,tmp_gwas_coloc_output_mqtl_df)
}

# joining coloc results to trait lists
joined_gwas_coloc_output_mqtl_df <- left_join(gwas_coloc_output_mqtl_df,cleaned_trait_category_names,by=c("V3"))
gwas_coloc_output_mqtl_df <- joined_gwas_coloc_output_mqtl_df

# identifying all mqtls that were included in the coloc analysis
included_mqtl <- gwas_coloc_output_mqtl_df %>% select(combination, cpg_id, variant_id) %>% unique()
nrow(included_mqtl)
length(unique(included_mqtl$variant_id))

# identifying mqtls that exhibited a coloc PP4>0.5 for at least one GWAS trait  
pp4_0.5_mqtl <- gwas_coloc_output_mqtl_df %>% filter(PP.H4.abf>0.5) %>% rename(Trait=V2) %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
unique_pp4_0.5_mqtl <- pp4_0.5_mqtl %>% select(combination, cpg_id, variant_id) %>% unique()
nrow(unique_pp4_0.5_mqtl)
length(unique(unique_pp4_0.5_mqtl$variant_id))

#######################################
ALL_JOINED_COLOC_I_MQTL <- inner_join(
  gwas_coloc_output_df %>% select(-PP.H0.abf,-PP.H1.abf,-PP.H2.abf,-PP.H3.abf) %>% rename(PP4_imqtl=PP.H4.abf),
  gwas_coloc_output_mqtl_df %>% select(-PP.H0.abf,-PP.H1.abf,-PP.H2.abf,-PP.H3.abf) %>% rename(PP4_mqtl=PP.H4.abf)
) %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))

# saving the list of mqtls that coloc with GWAS at imqtl
mqtl_coloc_at_imqtl <- pp4_0.5_mqtl$imqtl_id
save(mqtl_coloc_at_imqtl,file="/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/mqtl_coloc_at_imqtl.RData")

# identifying imQTLs uquiely identified by imQTL rather than mQTL
length(setdiff(pp4_0.5_imqtl$imqtl_id,pp4_0.5_mqtl$imqtl_id))
length(setdiff(pp4_0.5_mqtl$imqtl_id,pp4_0.5_imqtl$imqtl_id))


# joining gwas tables for coloc results
pp4_0.5_mqtl <- pp4_0.5_mqtl %>% select(-PP.H0.abf,-PP.H1.abf,-PP.H2.abf,-PP.H3.abf) %>% rename(PP4_mqtl=PP.H4.abf)
pp4_0.5_imqtl <- pp4_0.5_imqtl %>% select(-PP.H0.abf,-PP.H1.abf,-PP.H2.abf,-PP.H3.abf) %>% rename(PP4_imqtl=PP.H4.abf)
# basically joining for either significant as imQTL or mQTL
JOINED_COLOC_I_MQTL_pp4_0.5_InAtLeastOneQTL <- ALL_JOINED_COLOC_I_MQTL %>% filter(PP4_imqtl > 0.5 | PP4_mqtl > 0.5)

JOINED_COLOC_I_MQTL_pp4_0.5_InAtLeastOneQTL <- full_join(pp4_0.5_imqtl,pp4_0.5_mqtl)



#################################################################
# CREATING PLOT TO SEE HOW MANY PAIRS IMQTLS IMPROVE GWAS COLOC #
#################################################################
library(ggplot2)
library(dplyr)
library(tibble)

# Step 1: Define status categories
df <- JOINED_COLOC_I_MQTL_pp4_0.5_InAtLeastOneQTL %>%
  mutate(status = case_when(
    PP4_imqtl > PP4_mqtl ~ "PP4 of Interaction Term > Marginal Effect",
    PP4_imqtl < PP4_mqtl ~ "PP4 of Interaction Term < Marginal Effect",
    TRUE ~ "PP4 of Interaction Term = Marginal Effect"
  ))

# Step 2: Set up legend labels and colors
status_counts <- table(df$status)
desired_order <- c(
  "PP4 of Interaction Term > Marginal Effect",
  "PP4 of Interaction Term < Marginal Effect",
  "PP4 of Interaction Term = Marginal Effect"
)
present_statuses <- intersect(desired_order, names(status_counts))

label_map <- setNames(
  paste0(present_statuses, ": \n   ", as.integer(status_counts[present_statuses]), " colocalizing imQTL-GWAS pairs"),
  present_statuses
)

df <- df %>%
  mutate(
    legend_status = recode(status, !!!label_map),
    legend_status = factor(legend_status, levels = label_map)
  )

color_lookup <- c(
  "PP4 of Interaction Term > Marginal Effect" = "salmon",
  "PP4 of Interaction Term < Marginal Effect" = "skyblue",
  "PP4 of Interaction Term = Marginal Effect" = "gray40"
)
color_map <- setNames(color_lookup[present_statuses], label_map[present_statuses])

# Step 3: Assign quadrants using > 0.5 and <= 0.5
df <- df %>%
  mutate(quadrant = case_when(
    PP4_mqtl > 0.5 & PP4_imqtl > 0.5 ~ "top_right",
    PP4_mqtl <= 0.5 & PP4_imqtl > 0.5 ~ "top_left",
    PP4_mqtl > 0.5 & PP4_imqtl <= 0.5 ~ "bottom_right",
    PP4_mqtl <= 0.5 & PP4_imqtl <= 0.5 ~ "bottom_left",
    TRUE ~ NA_character_
  ))

# Step 4: Count quadrants and fill missing ones
quad_counts <- df %>%
  count(quadrant) %>%
  tibble::deframe()

all_quads <- c("top_left", "top_right", "bottom_right", "bottom_left")
quad_counts <- setNames(quad_counts[all_quads], all_quads)
quad_counts[is.na(quad_counts)] <- 0

# Step 5: Plot with updated quadrant labels
output_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/GWAS_coloc_scatterplot_pp4_imqtl_mqtl.png"

png(output_path, width = 5.5, height = 5.5, units = "in", res = 300)
ggplot(df, aes(x = PP4_mqtl, y = PP4_imqtl, color = legend_status)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "black") +
  geom_vline(xintercept = 0.5, color = "black") +
  annotate("text", x = 0.25, y = 0.75,
           label = paste0("bold('", quad_counts["top_left"], " pairs')"),
           size = 4.5, parse = TRUE) +
  annotate("text", x = 0.75, y = 0.75,
           label = paste0("bold('", quad_counts["top_right"], " pairs')"),
           size = 4.5, parse = TRUE) +
  annotate("text", x = 0.75, y = 0.25,
           label = paste0("bold('", quad_counts["bottom_right"], " pairs')"),
           size = 4.5, parse = TRUE) +
  labs(
    x = "PP4 (Coloc using marginal effects)",
    y = "PP4 (Coloc using interaction terms)",
    title = "",
    color = ""
  ) +
  scale_color_manual(values = color_map) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
dev.off()



##################################################################
# PIE CHART: INTERACTION VS MARGINAL COLOCALIZATION CATEGORIES  #
##################################################################

library(dplyr)
library(ggplot2)

# Input data
df <- JOINED_COLOC_I_MQTL_pp4_0.5_InAtLeastOneQTL

# Summarize and categorize
pp4_summary <- df %>%
  group_by(imqtl_id) %>%
  summarise(
    max_PP4_interaction = max(PP4_imqtl, na.rm = TRUE),
    max_PP4_marginal = max(PP4_mqtl, na.rm = TRUE)
  ) %>%
  mutate(
    category = case_when(
      max_PP4_interaction > 0.5 & max_PP4_marginal <= 0.5 ~ "Interaction-only colocalization",
      max_PP4_interaction <= 0.5 & max_PP4_marginal > 0.5 ~ "Marginal-only colocalization",
      max_PP4_interaction > max_PP4_marginal ~ "Both colocalized (higher in Interaction)",
      max_PP4_interaction < max_PP4_marginal ~ "Both colocalized (higher in Marginal)",
      TRUE ~ "Equal PP4"
    )
  )

# Set category order and final color palette
category_colors <- c(
  "Interaction-only colocalization" = "salmon",        # darker red
  "Marginal-only colocalization" = "#007FFF",          # darker blue
  "Both colocalized (higher in Interaction)" = "#FFD1DC",  # lighter red
  "Both colocalized (higher in Marginal)" = "lightblue",     # lighter blue
  "Equal PP4" = "gray50"
)

pp4_summary$category <- factor(pp4_summary$category, levels = c(
  "Interaction-only colocalization",
  "Both colocalized (higher in Interaction)",
  "Equal PP4",
  "Both colocalized (higher in Marginal)",
  "Marginal-only colocalization"
))

# Count frequencies
category_counts <- pp4_summary %>%
  count(category) %>%
  mutate(
    percent = round(100 * n / sum(n), 1),
    label = paste0(n, " (", percent, "%)")
  )

# Create pie chart
pie_plot <- ggplot(category_counts, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  geom_text(
    aes(x = 2, label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4.5
  ) +
  scale_fill_manual(values = category_colors) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# Save high-resolution PNG
ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/GWAS_coloc_pie_pp4_interaction_vs_marginal.png",
  plot = pie_plot,
  width = 7,
  height = 7,
  units = "in",
  dpi = 300
)

