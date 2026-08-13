library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

##################################
# importing eqtl colocalizations #
##################################
all_eqtl_imqtl_coloc <- fread("/gpfs/data/pierce-lab/james.li/imQTL/output/EQTL_coloc/all_eqtl_imqtl_coloc.tsv") %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
eqtl_coloc_list <- unique((all_eqtl_imqtl_coloc %>% filter(coloc_0.5_all==TRUE))$imqtl_id)

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
nrow(included_imqtl)
length(unique(included_imqtl$variant_id))

# identifying imqtls that exhibited a coloc PP4>0.5 for at least one GWAS trait  
pp4_0.5_imqtl <- gwas_coloc_output_df %>% filter(PP.H4.abf>0.5) %>% rename(Trait=V2) 
unique_pp4_0.5_imqtl <- pp4_0.5_imqtl %>% select(combination, cpg_id, variant_id) %>% unique()
nrow(unique_pp4_0.5_imqtl)
length(unique(unique_pp4_0.5_imqtl$variant_id))

#################################
# creating a unique combination, cpg_id, variant_id list for imQTLs that colocalized with GWAS signals
pp4_0.5_imqtl <- pp4_0.5_imqtl %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
gwas_coloc_list <- unique((pp4_0.5_imqtl)$imqtl_id)

# identifying imqtls with multiple colocalization signals
counts <- data.frame(pp4_0.5_imqtl %>% select(imqtl_id) %>% table()) %>% arrange(desc(Freq)) 

############################
# plot histogram for number of colocalizing signals
############################

# Set output file path
output_path <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/GWAS_coloc_multiple_hit_histogram.png"

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
pp4_counts <- pp4_0.5_imqtl %>%
  count(Category, Trait) %>%
  mutate(
    Trait = factor(Trait),
    Category = factor(Category)
  )

# Step 2: Plot with facet_grid to simulate grouped spacing
png("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/GWAS_coloc_hist.png", 
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


######################################
# Plot stacked histogram (by tissue) #
######################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Create tissue lookup vector
tissue_map <- c(
  "colon"    = "Colon",
  "lung"     = "Lung",
  "prostate" = "Prostate",
  "wb"       = "Whole Blood",
  "ovary"    = "Ovary"
)

# Custom tissue colors
tissue_colors <- c(
  "Lung"        = "yellowgreen",
  "Colon"       = "sienna",
  "Ovary"       = "pink3",
  "Prostate"    = "lightgray",
  "Whole Blood" = "magenta3"
)

# Count colocalizations by category, trait, and tissue
pp4_counts <- pp4_0.5_imqtl %>%
  mutate(
    Category = ifelse(
      Category == "Cardiometabolic",
      "Cardiometabolic",
      Category
    ),
    Category = ifelse(
      Category == "Endocrine system disease",
      "Endo.",
      Category
    ),
    Category = ifelse(
      Category == "Skeletal system disease",
      "Skeletal",
      Category
    ),
    Category = ifelse(
      Category == "Digestive system disease",
      "Digest.  ",
      Category
    ),
    Category = ifelse(
      Category == "Aging",
      "Aging ",
      Category
    ),
    Tissue = recode(tissue, !!!tissue_map)
  ) %>%
  count(Category, Trait, Tissue) %>%
  mutate(
    Trait = factor(Trait),
    Category = factor(Category),
    Tissue = factor(
      Tissue,
      levels = names(tissue_colors)
    )
  ) %>%
  complete(
    nesting(Category, Trait),
    Tissue = factor(
      names(tissue_colors),
      levels = names(tissue_colors)
    ),
    fill = list(n = 0)
  )

png(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/GWAS_coloc_stacked_hist.png",
  res = 300,
  units = "in",
  height = 4.8,
  width = 13
)

ggplot(
  pp4_counts,
  aes(
    x = Trait,
    y = n,
    fill = Tissue
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 0.35
  ) +
  facet_grid(
    . ~ Category,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    x = NULL,
    y = "Number of colocalizations"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    
    # Horizontal category labels above each facet
    strip.background = element_blank(),
    strip.clip = "off",
    strip.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 8.5,
      face = "bold",
      margin = margin(
        t = 4,
        r = 3,
        b = 8,
        l = 3
      )
    ),
    
    # More space between category groups
    panel.spacing.x = unit(
      0.9,
      "lines"
    ),
    
    # Remove panel boxes
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Light horizontal reference grid
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.3
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    plot.margin = margin(
      t = 12,
      r = 12,
      b = 10,
      l = 30
    ),
    
    legend.text = element_text(
      size = 7.5
    ),
    legend.title = element_blank()
  ) +
  scale_fill_manual(
    values = tissue_colors,
    breaks = names(tissue_colors),
    limits = names(tissue_colors),
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )
dev.off()




######################
# Plot arc pie chart #
######################
library(ggplot2)
library(dplyr)

# Replace with real vectors
length(eqtl_coloc_list)
length(gwas_coloc_list)
gwas_coloc_output_df <- gwas_coloc_output_df %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep="."))
all_tested_imqtl_list <- unique(gwas_coloc_output_df$imqtl_id)
length(all_tested_imqtl_list)

# Ensure uniqueness
all <- unique(all_tested_imqtl_list)
gwas <- unique(gwas_coloc_list)
eqtl <- unique(eqtl_coloc_list)

# Define groups
both_coloc <- intersect(gwas, eqtl)
gwas_only <- setdiff(gwas, eqtl)
eqtl_only <- setdiff(eqtl, gwas)
non_coloc <- setdiff(all, union(gwas, eqtl))

# Print out proportion
length(both_coloc)
length(gwas_only)
length(eqtl_only)
length(non_coloc)

length(both_coloc)/length(all_tested_imqtl_list)
length(gwas_only)/length(all_tested_imqtl_list)
length(eqtl_only)/length(all_tested_imqtl_list)
length(non_coloc)/length(all_tested_imqtl_list)

# Define correct angular order
arc_df <- data.frame(
  category = c(
    "Non-colocalizing imQTL",
    "Colocalize with GWAS only",
    "Colocalize with GWAS+eQTL",
    "Colocalize with eQTL only"
  ),
  count = c(
    length(non_coloc),
    length(gwas_only),
    length(both_coloc),
    length(eqtl_only)
  ),
  r0 = c(0.45, 0.5, 0.55, 0.5)-0.3,
  r = c(0.5, 0.55, 0.6, 0.55)-0.3
)

# Compute arc angles
total <- sum(arc_df$count)
arc_df <- arc_df %>%
  mutate(
    frac = count / total,
    start = cumsum(c(0, head(frac, -1))) * 2 * pi,
    end = cumsum(frac) * 2 * pi,
    mid = (start + end) / 2,
    label_x = 1.1 * cos(mid),
    label_y = 1.1 * sin(mid),
    label = paste0(category, "\n(n = ", count, ")")
  )

# Create arc polygons
create_arc_polygon <- function(start_angle, end_angle, r0, r, n = 100, category) {
  theta_outer <- seq(start_angle, end_angle, length.out = n)
  theta_inner <- rev(theta_outer)
  x_outer <- r * cos(theta_outer)
  y_outer <- r * sin(theta_outer)
  x_inner <- r0 * cos(theta_inner)
  y_inner <- r0 * sin(theta_inner)
  data.frame(
    x = c(x_outer, x_inner),
    y = c(y_outer, y_inner),
    category = category
  )
}

# Generate all arc polygons
poly_data <- do.call(rbind, lapply(1:nrow(arc_df), function(i) {
  create_arc_polygon(
    start_angle = arc_df$start[i],
    end_angle = arc_df$end[i],
    r0 = arc_df$r0[i],
    r = arc_df$r[i],
    category = arc_df$category[i]
  )
}))

# Plot
p <- ggplot() +
  geom_polygon(data = poly_data,
               aes(x, y, fill = category, group = category),
               color = "black") +
  geom_text(data = arc_df,
            aes(x = label_x, y = label_y, label = label,
                hjust = ifelse(cos(mid) > 0, 0, 1)),
            size = 3.5) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    text = element_text(size = 10)
  ) +
  scale_fill_manual(values = c(
    "Colocalize with GWAS+eQTL" = "#cba6f7",  # light purple
    "Colocalize with GWAS only" = "#87ceeb",  # sky blue
    "Colocalize with eQTL only" = "#fa8072",  # salmon
    "Non-colocalizing imQTL"    = "#999999"   # gray
  ))

# Save
ggsave(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/arc_chart.png",
  plot = p, width = 8, height = 8, dpi = 1200
)

###########################################
# Plot stacked arc pie chart (by tissue)  #
###########################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Replace with real vectors
length(eqtl_coloc_list)
length(gwas_coloc_list)

gwas_coloc_output_df <- gwas_coloc_output_df %>%
  mutate(
    imqtl_id = paste(combination, cpg_id, variant_id, sep = ".")
  )

all_tested_imqtl_list <- unique(gwas_coloc_output_df$imqtl_id)
length(all_tested_imqtl_list)

# Ensure uniqueness
all <- unique(all_tested_imqtl_list)
gwas <- unique(gwas_coloc_list)
eqtl <- unique(eqtl_coloc_list)

# Define groups
both_coloc <- intersect(gwas, eqtl)
gwas_only <- setdiff(gwas, eqtl)
eqtl_only <- setdiff(eqtl, gwas)
non_coloc <- setdiff(all, union(gwas, eqtl))

# Print out proportions
length(both_coloc)
length(gwas_only)
length(eqtl_only)
length(non_coloc)

length(both_coloc) / length(all_tested_imqtl_list)
length(gwas_only) / length(all_tested_imqtl_list)
length(eqtl_only) / length(all_tested_imqtl_list)
length(non_coloc) / length(all_tested_imqtl_list)

# Tissue lookup
tissue_map <- c(
  "colon"    = "Colon",
  "lung"     = "Lung",
  "prostate" = "Prostate",
  "wb"       = "Whole Blood",
  "ovary"    = "Ovary"
)

# Custom tissue colors
tissue_colors <- c(
  "Lung"        = "yellowgreen",
  "Colon"       = "sienna",
  "Ovary"       = "pink3",
  "Prostate"    = "lightgray",
  "Whole Blood" = "magenta3"
)

# Obtain one tissue assignment for each tested imQTL
imqtl_tissue_df <- gwas_coloc_output_df %>%
  distinct(imqtl_id, tissue) %>%
  mutate(
    Tissue = recode(tissue, !!!tissue_map),
    Tissue = factor(
      Tissue,
      levels = names(tissue_colors)
    )
  )

# Assign each imQTL to a colocalization category
arc_long_df <- imqtl_tissue_df %>%
  mutate(
    category = case_when(
      imqtl_id %in% non_coloc  ~ "Non-colocalizing imQTL",
      imqtl_id %in% gwas_only  ~ "Colocalize with GWAS only",
      imqtl_id %in% both_coloc ~ "Colocalize with GWAS+eQTL",
      imqtl_id %in% eqtl_only  ~ "Colocalize with eQTL only"
    )
  ) %>%
  filter(!is.na(category)) %>%
  count(category, Tissue, name = "count") %>%
  complete(
    category,
    Tissue = factor(
      names(tissue_colors),
      levels = names(tissue_colors)
    ),
    fill = list(count = 0)
  )

# Define correct angular category order
category_order <- c(
  "Non-colocalizing imQTL",
  "Colocalize with GWAS only",
  "Colocalize with GWAS+eQTL",
  "Colocalize with eQTL only"
)

arc_long_df <- arc_long_df %>%
  mutate(
    category = factor(category, levels = category_order)
  ) %>%
  arrange(category, Tissue)

# Define radius for each category
radius_df <- data.frame(
  category = factor(category_order, levels = category_order),
  r0 = c(0.45, 0.50, 0.55, 0.50) - 0.30,
  r  = c(0.50, 0.55, 0.60, 0.55) - 0.30
)

arc_long_df <- arc_long_df %>%
  left_join(radius_df, by = "category")

# Compute arc angles for tissue-specific segments
total <- sum(arc_long_df$count)

arc_long_df <- arc_long_df %>%
  mutate(
    frac = count / total,
    start = cumsum(lag(frac, default = 0)) * 2 * pi,
    end = cumsum(frac) * 2 * pi
  )

# Compute one label per colocalization category
arc_labels <- arc_long_df %>%
  group_by(category) %>%
  summarise(
    count = sum(count),
    start = min(start),
    end = max(end),
    .groups = "drop"
  ) %>%
  mutate(
    mid = (start + end) / 2,
    label_x = 1.1 * cos(mid),
    label_y = 1.1 * sin(mid),
    label = paste0(category, "\n(n = ", count, ")")
  )

# Create arc polygons
create_arc_polygon <- function(
    start_angle,
    end_angle,
    r0,
    r,
    n = 100,
    category,
    Tissue,
    segment_id
) {
  theta_outer <- seq(start_angle, end_angle, length.out = n)
  theta_inner <- rev(theta_outer)
  
  x_outer <- r * cos(theta_outer)
  y_outer <- r * sin(theta_outer)
  x_inner <- r0 * cos(theta_inner)
  y_inner <- r0 * sin(theta_inner)
  
  data.frame(
    x = c(x_outer, x_inner),
    y = c(y_outer, y_inner),
    category = category,
    Tissue = Tissue,
    segment_id = segment_id
  )
}

# Generate tissue-specific arc polygons
poly_data <- do.call(
  rbind,
  lapply(seq_len(nrow(arc_long_df)), function(i) {
    
    # Skip zero-count tissue segments
    if (arc_long_df$count[i] == 0) {
      return(NULL)
    }
    
    create_arc_polygon(
      start_angle = arc_long_df$start[i],
      end_angle = arc_long_df$end[i],
      r0 = arc_long_df$r0[i],
      r = arc_long_df$r[i],
      category = as.character(arc_long_df$category[i]),
      Tissue = as.character(arc_long_df$Tissue[i]),
      segment_id = i
    )
  })
)

# Plot
p <- ggplot() +
  geom_polygon(
    data = poly_data,
    aes(
      x = x,
      y = y,
      fill = Tissue,
      group = segment_id
    ),
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = arc_labels,
    aes(
      x = label_x,
      y = label_y,
      label = label,
      hjust = ifelse(cos(mid) > 0, 0, 1)
    ),
    size = 3.5
  ) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    text = element_text(size = 10)
  ) +
  scale_fill_manual(
    values = tissue_colors,
    breaks = names(tissue_colors),
    limits = names(tissue_colors),
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

# Save
ggsave(
  "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/arc_chart_stacked.png",
  plot = p,
  width = 8,
  height = 8,
  dpi = 1200
)

###################################
# simple pie chart of coloc vs not coloc
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
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/coloc_piechart_gwas_imqtl_default_prior.png",
  plot = p,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)




###########################################################
# importing imQTLs with colocalized mQTL marginal effects #
###########################################################
# load("/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/mqtl_coloc_at_imqtl.RData")




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
OUTPUT_TABLE_GWAS_COLOC <- pp4_0.5_imqtl %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep=".")) %>%
  separate(combination, into = c("tissue", "celltype"), sep = "_") %>%
  mutate(
    tissue = tissue_map[tissue],
    celltype = celltype_map[celltype]
  ) 
# writing out this table
write.table(OUTPUT_TABLE_GWAS_COLOC,file="/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/OUTPUT_TABLE_GWAS_COLOC.tsv",quote=F,sep="\t",row.names=F,col.names=T)


#################################
# WRITING OUT ALL COLOC RESULTS #
#################################
# Apply the mapping
ALL_OUTPUT_TABLE_GWAS_COLOC <- gwas_coloc_output_df %>% rename(Trait=V2) %>% mutate(imqtl_id=paste(combination,cpg_id,variant_id,sep=".")) %>%
  separate(combination, into = c("tissue", "celltype"), sep = "_") %>%
  mutate(
    tissue = tissue_map[tissue],
    celltype = celltype_map[celltype]
  ) 
# writing out this table
write.table(ALL_OUTPUT_TABLE_GWAS_COLOC,file="/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/ALL_OUTPUT_TABLE_GWAS_COLOC.tsv",quote=F,sep="\t",row.names=F,col.names=T)










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
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/revision/PP_HIST.png",
  plot = pp_hist_plot,
  width = 6, height = 8, dpi = 300
)
