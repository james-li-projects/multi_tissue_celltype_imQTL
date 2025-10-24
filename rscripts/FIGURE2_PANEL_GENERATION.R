library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(ggnewscale) 
library(dplyr)
library(tidyr)
library(purrr)

##############################
# set working directory
setwd("/gpfs/data/pierce-lab/james.li/imQTL/output/regional_plots")

##############################
# User-supplied values [class 1]
query_phenotype <- "cg24996280"
query_variant <- "chr8:142399832:G:A"
query_rsid <- "rs68095309"
dataset="GTEx"
current_tissue="colon"
tmp_combination="colon_Epi"

##############################
# User-supplied values [class 2]
query_phenotype <- "cg06643156"
query_variant <- "chr3:58493627:A:G"
query_rsid <- "rs73087706"
dataset="GTEx"
current_tissue="colon"
tmp_combination="colon_Stromal"

##############################
# User-supplied values [class 3]
query_phenotype <- "cg02044966"
query_variant <- "chr17:65907762:T:G"
query_rsid <- "rs9906405"
dataset="HEALS"
current_tissue="wb"
tmp_combination="wb_Neutro"


##############################
# User-supplied values -- interaction coloc #1
query_phenotype <- "cg04093349"
query_variant <- "chr15:66716105:A:G"
query_rsid <- "rs56173559"
dataset="GTEx"
current_tissue="colon"
tmp_combination="colon_Epi"

##############################
# User-supplied values -- interaction coloc #2
query_phenotype <- "cg13848087"
query_variant <- "chr12:89456163:A:G"
query_rsid <- "rs10858866"
dataset="GTEx"
current_tissue="colon"
tmp_combination="colon_Stromal"


##############################
# User-supplied values -- interaction coloc #3
query_phenotype <- "cg20356878"
query_variant <- "chr3:121995837:G:A"
query_rsid <- "rs9816720"
dataset="GTEx"
current_tissue="lung"
tmp_combination="lung_Mono"

##############################
# User-supplied values -- interaction coloc #4
query_phenotype <- "cg19255693"
query_variant <- "chr20:48927097:A:C"
query_rsid <- "rs35393280"
dataset="GTEx"
current_tissue="lung"
tmp_combination="lung_Macro"


# manual querying of results for examples
load("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.RData")
wide_parsed_imqtl %>% filter(tissue=="wb") %>% filter(P_MarginalEffect>0.05,P_InteractionEffect<1e-11) %>% arrange(desc(P_MarginalEffect),P_InteractionEffect) %>% head(10)
wide_parsed_imqtl %>% filter(phenotype_id=="cg02044966")


##############################
# importing pvar
pvar<-fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/data/",dataset,"/genetic_data/processed_genetic_data_chrprefix.pvar"))

# importing tissue DNA methylation bed file
input_dnam_bed <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",dataset,"/",current_tissue,".bed"),sep="\t") 
# importing tissue covariates files
cov_df <- data.frame(t(read.table(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",dataset,"/processed_covariates_",current_tissue,".txt"),sep="\t"))) 
if (dataset=="GTEx") {
  cov_df <- cov_df %>% mutate(Sample_Name=gsub("\\.","-",rownames(cov_df)))
} else if (dataset=="HEALS") {
  cov_df <- cov_df %>% mutate(Sample_Name=gsub("^X","",rownames(cov_df)))
}
# write out tissue list of samples
write.table(data.frame(0,cov_df$Sample_Name),file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list",quote=F,row.names=F,col.names=F,sep="\t")
    
# importing CT interaction file
ct_df <- read.table(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",dataset,"/processed_interactions_",tmp_combination,".txt"),sep="\t") %>% rename(Sample_Name=V1,CT=V2)

# importing imQTL results
tmp_imqtl_df <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/",dataset,"/imQTL/top_assoc/tensorQTL_imQTL_",tmp_combination,".cis_qtl_top_assoc.txt.gz")) %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni")) %>% filter(pval_adj_bonf<0.05)

# Get the row number
b <- which(tmp_imqtl_df$phenotype_id == query_phenotype & tmp_imqtl_df$variant_id == query_variant)

# initializing target cpg and lead variant
target_cpg = tmp_imqtl_df$phenotype_id[b]
lead_variant_raw = tmp_imqtl_df$variant_id[b]

###############################
# reformatting lead variant with periods for later matching
lead_variant_period<-gsub("\\:",".",lead_variant_raw)

# importing all variants that should be utilized in regional plot
num_chr<-as.numeric(gsub("chr","",(stringr::str_split(lead_variant_raw,pattern="\\:"))[[1]][1]))
num_pos<-as.numeric(gsub("chr","",(stringr::str_split(lead_variant_raw,pattern="\\:"))[[1]][2]))
window_variants <- pvar %>% filter(`#CHROM`==paste0("chr",num_chr)) %>% filter(POS>(num_pos-250000),POS<(num_pos+250000))
write.table(window_variants %>% select(ID),file="/gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list",quote=F,row.names=F,col.names=F,sep="\t")

# filering DNAm bed file 
dnam_bed <- input_dnam_bed %>% filter(phenotype_id==target_cpg) %>% head(1)
t_dnam_bed <- data.frame(t(dnam_bed))
t_dnam_bed <- t_dnam_bed %>% slice(-c(1:4)) 
cpg_df <- t_dnam_bed %>% rename(mval=t.dnam_bed.) %>% mutate(Sample_Name=rownames(t_dnam_bed)) %>% select(Sample_Name,everything())

# importing variant data
system(paste0("module load plink/2.0; plink2 -pfile /gpfs/data/pierce-lab/james.li/imQTL/data/",dataset,"/genetic_data/processed_genetic_data_chrprefix --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list --keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list --maf 0.05 --export Av --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window"))
traw <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_window.traw") %>% separate(SNP,into=c("T1","T2","T3","T4"),remove=F) %>% mutate(across(11:ncol(.), ~ ifelse(COUNTED != T4, 2 - ., .))) %>% select(-c(paste0(rep("T","4"),1:4)))

# obtaining df genetic variants
variant_df <- 
  data.frame(t(traw))
newcolnames <- variant_df["SNP",]
variant_df <- variant_df %>%
  rename_with(~ as.character( newcolnames ), everything())
# Filter out the rows to be removed
rows_to_remove <- c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT")
variant_df <- variant_df %>%
  filter(!rownames(.) %in% rows_to_remove)
# adding sample name column
variant_df <- variant_df %>% mutate(Sample_Name = gsub("^0_","",rownames(.))) %>% select(Sample_Name,everything())

########################
# Join the data frames by "Sample_Name" using inner_join
reg_df <- cpg_df %>%
  inner_join(ct_df, by = "Sample_Name") %>%
  inner_join(cov_df, by = "Sample_Name") %>%
  inner_join(variant_df, by = "Sample_Name") %>% select(-Sample_Name)
colnames(reg_df) <- gsub("\\:",".",colnames(reg_df))

# Convert all predictors to numeric
convert_to_numeric <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      return(col)
    } else {
      # Convert non-numeric columns to numeric
      return(as.numeric(col))
    }
  })
  return(df)
}
# Ensure all predictors are numeric
reg_df <- convert_to_numeric(reg_df)



################################
# Generate interaction boxplot #
################################

# Extract unprocessed names from tmp_combination
unprocessed_tissue <- str_split(tmp_combination, pattern = "_")[[1]][1]
unprocessed_celltype <- str_split(tmp_combination, pattern = "_")[[1]][2]

# Create lookup vectors
tissue_map <- c(
  "breast" = "Breast", "colon" = "Colon", "lung" = "Lung", "kidney" = "Kidney",
  "prostate" = "Prostate", "wb" = "Whole Blood", "ovary" = "Ovary"
)
celltype_map <- c(
  "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
  "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "SM" = "Smooth muscle cell",
  "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
  "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
  "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
  "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell",
  "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
)

# Map to pretty labels
processed_tissue <- tissue_map[[unprocessed_tissue]]
processed_celltype <- celltype_map[[unprocessed_celltype]]
processed_tmp_combination <- paste0(processed_tissue, " - ", processed_celltype)

# Prepare data for interaction boxplot
interaction_boxplot_df <- reg_df %>%
  select(mval, CT, all_of(lead_variant_period)) %>%
  rename_with(~ "dosage", .cols = 3) %>%
  mutate(dosage = round(dosage, digits = 0)) %>%
  filter(dosage %in% c(0, 1, 2))

# Format genotype labels
alleles <- unlist(str_split(lead_variant_period, pattern = "\\."))[3:4]
hom_rec_str <- paste0(alleles[1], "/", alleles[1])
het_str     <- paste0(alleles[1], "/", alleles[2])
hom_dom_str <- paste0(alleles[2], "/", alleles[2])

interaction_boxplot_df <- interaction_boxplot_df %>%
  mutate(Genotype = case_when(
    dosage == 0 ~ hom_rec_str,
    dosage == 1 ~ het_str,
    dosage == 2 ~ hom_dom_str
  )) %>%
  filter(!is.na(Genotype)) %>% mutate(Genotype=factor(Genotype,levels=c(hom_rec_str,het_str,hom_dom_str)))

# Cell type proportion groups
interaction_boxplot_df <- interaction_boxplot_df %>%
  mutate(`Cell Type Proportion` = ifelse(
    CT > median(CT),
    paste0("Upper 50% (", processed_celltype, ")"),
    paste0("Lower 50% (", processed_celltype, ")")
  ))

# Add "All Individuals" group
copy_df <- interaction_boxplot_df %>%
  mutate(`Cell Type Proportion` = "All Individuals")

interaction_boxplot_df <- bind_rows(interaction_boxplot_df, copy_df)

# Convert genotype to numeric for regression line
interaction_boxplot_df <- interaction_boxplot_df %>%
  mutate(Genotype_numeric = as.numeric(as.factor(Genotype)))

################################
# start of plotting code
library(dplyr)
library(ggplot2)
library(broom)

x_label_pos <- max(interaction_boxplot_df$Genotype_numeric, na.rm = TRUE) - 0.6
y_label_pos <- min(interaction_boxplot_df$mval, na.rm = TRUE) + 0.1

# Initialize an empty list to store results
label_list <- list()

# Get unique cell type proportions
cell_type_levels <- unique(interaction_boxplot_df$`Cell Type Proportion`)

# Loop through each cell type proportion
for (i in seq_along(cell_type_levels)) {
  current_level <- cell_type_levels[i]
  
  # Subset data
  subset_df <- interaction_boxplot_df[interaction_boxplot_df$`Cell Type Proportion` == current_level, ]
  
  # Fit regression
  model <- lm(mval ~ dosage, data = subset_df)
  tidy_model <- broom::tidy(model)
  
  # Extract stats
  beta <- tidy_model$estimate[2]
  pval <- tidy_model$p.value[2]
  
  # Format labels
  beta_label <- sprintf("\u03B2 = %.2f", beta)
  pval_label <- if (pval < 0.05) {
    paste0("p = ", formatC(pval, format = "e", digits = 1))
  } else {
    paste0("p = ", sprintf("%.2f", pval))
  }
  
  # Store result
  label_list[[i]] <- data.frame(
    `Cell Type Proportion` = current_level,
    label = paste(beta_label, pval_label, sep = "\n"),
    x = x_label_pos,
    y = y_label_pos
  )
}

# Combine all results into one data frame
regression_labels <- do.call(rbind, label_list)

# Reattach facet variable
regression_labels$`Cell Type Proportion` <- unique(interaction_boxplot_df$`Cell Type Proportion`)

interact_boxplot_p <- ggplot(interaction_boxplot_df, aes(x = Genotype, y = mval)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = Genotype)) +
  geom_smooth(
    aes(x = Genotype_numeric, y = mval),
    method = "lm", se = FALSE, color = "darkgray", linewidth = 1
  ) +
  geom_text(
    data = regression_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0,
    size = 4
  ) +
  facet_wrap(~ `Cell Type Proportion`, nrow = 1) +
  theme_classic() +
  labs(
    x = "Genotype",
    y = "DNAm (INT M-value)",
    title = paste("DNAm at", target_cpg, "by", query_rsid, "and", tolower(processed_celltype), "proportion in", tolower( processed_tissue))
  ) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  )

ggsave(
  filename = paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/regional_plots/figure2/boxplot_",
                    tmp_combination, "_", lead_variant_period, "_", target_cpg, ".png"),
  plot = interact_boxplot_p,
  width = 9.5,
  height = 4
)







############################
# GENERATING REGIONAL PLOT #
############################
# Identify column indices of reg_df that start with "chr"
chr_indices <- grep("^chr", colnames(reg_df))

# Identify the index of the 'lead_variant'
lead_variant_index <- which(gsub("\\:",".",colnames(reg_df)) == lead_variant_period)

# Get indices of columns starting with "chr" but not equal to 'lead_variant_raw'
evaluated_variant_indices <- setdiff(chr_indices, lead_variant_index)

#######################################################
# Function to run regression and extract coefficients #
#######################################################
run_regression <- function(df, outcome, predictors) {
  # Construct the formula
  formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  
  # Generate model matrix to check rank deficiency
  model_matrix <- model.matrix(formula, data = df)
  model_rank <- qr(model_matrix)$rank  # Get rank of matrix
  expected_rank <- length(predictors) + 1  # Expected rank (including intercept)
  
  # If rank is lower than expected, return NA for all coefficients
  if (model_rank < expected_rank) {
    coef_names <- c("(Intercept)", predictors)  # Ensure correct variable names
    return(matrix(NA, nrow = length(coef_names), ncol = 4, 
                  dimnames = list(coef_names, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))))
  }
  
  # Run linear regression
  model <- lm(formula, data = df)
  
  # Extract coefficients from summary if no collinearity issue
  return(summary(model)$coefficients)
}

################################################
# INITIALIZING DATA FRAME TO STORE ALL RESULTS #
################################################
results <- data.frame()

###################################################
# OBTAINING MARGINAL, MAIN, AND INTERACTION TERMS #
###################################################
# initializing list of baseline predictors and the outcome variable
tmp_outcome <- "mval"
tmp_baseline_predictors <- setdiff(colnames(reg_df[-c(chr_indices)]),"mval")
for (current_variant_col_index in chr_indices) {
  # obtaining name of variant being analyzed
  current_variant_col_name <- colnames(reg_df)[current_variant_col_index]
  
  ############
  # Marginal #
  ############
  # obtaining list of predictors
  tmp_predictors <- setdiff(c(current_variant_col_name,tmp_baseline_predictors),"CT")
  # obtaining regression value for current variant 
  regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
  tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Marginal Effect")
  results<-rbind(results,tmp_results)
  
  ################################
  # Interaction and Main Effects #
  ################################
  # obtaining list of predictors
  tmp_interaction_term <- paste0(current_variant_col_name,"*CT")
  tmp_predictors <- c(current_variant_col_name,tmp_baseline_predictors,tmp_interaction_term)
  # obtaining regression value for current variant 
  regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
  
  ###############
  # Interaction #
  ###############
  tmp_results <- cbind(regression_results[gsub("\\*",":",tmp_interaction_term),],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Interaction Effect")
  results<-rbind(results,tmp_results)
  
  ########
  # Main #
  ########
  tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Main Effect")
  results<-rbind(results,tmp_results)
}

###############################################################
# OBTAINING CONDITIONAL MARGINAL, MAIN, AND INTERACTION TERMS #
###############################################################
# obtaining marginal effects
tmp_outcome <- "mval"
tmp_baseline_predictors <- setdiff(colnames(reg_df[-c(chr_indices)]),"mval")
for (current_variant_col_index in evaluated_variant_indices) {
  # obtaining name of variant being analyzed
  current_variant_col_name <- colnames(reg_df)[current_variant_col_index]
  
  ############
  # Marginal #
  ############
  # obtaining list of predictors including the lead variant
  tmp_predictors <- c(setdiff(c(current_variant_col_name,tmp_baseline_predictors),"CT"),lead_variant_period)
  # obtaining regression value for current variant 
  regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
  tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Marginal Effect (conditional)")
  results<-rbind(results,tmp_results)
  
  ################################
  # Interaction and Main Effects #
  ################################
  # obtaining list of predictors including interaction terms with both the lead and target variants
  tmp_interaction_term_current_variant <- paste0(current_variant_col_name,"*CT")
  tmp_interaction_term_lead_variant <- paste0(lead_variant_period,"*CT")
  tmp_predictors <- c(current_variant_col_name,lead_variant_period,tmp_baseline_predictors,tmp_interaction_term_current_variant,tmp_interaction_term_lead_variant)
  # obtaining regression value for current variant 
  regression_results <- data.frame(run_regression(df = reg_df,outcome = tmp_outcome,predictors = tmp_predictors))
  
  ###############
  # Interaction #
  ###############
  tmp_results <- cbind(regression_results[gsub("\\*",":",tmp_interaction_term_current_variant),],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Interaction Effect (conditional)")
  results<-rbind(results,tmp_results)
  
  ########
  # Main #
  ########
  tmp_results <- cbind(regression_results[current_variant_col_name,],variant_id=current_variant_col_name,phenotype_id=target_cpg,class="Main Effect (conditional)")
  results<-rbind(results,tmp_results)
}


####################################
# Add r2 correlations between SNPs #
####################################
# computing LD
max_snp_window<-nrow(window_variants)
system(paste0("module load plink/1.9; plink -bfile /gpfs/data/pierce-lab/james.li/imQTL/data/",dataset,"/genetic_data/processed_genetic_data_chrprefix_bfile --extract /gpfs/data/pierce-lab/james.li/imQTL/tmp/variant.list --keep /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_samp.list --maf 0.05 --r2 --ld-window-r2 0 --ld-window ",max_snp_window," --out /gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_LD"))

# import pairwise correlations
pairwise_cor <- fread("/gpfs/data/pierce-lab/james.li/imQTL/tmp/tmp_LD.ld")

# Filter the rows where either SNP_A or SNP_B equals lead_variant_raw
filtered_df <- pairwise_cor[pairwise_cor$SNP_A == lead_variant_raw | pairwise_cor$SNP_B == lead_variant_raw, ]

# Create a new column 'OTHER_VARIANT'
filtered_df$OTHER_VARIANT <- ifelse(filtered_df$SNP_A == lead_variant_raw, filtered_df$SNP_B, filtered_df$SNP_A)

# finalizing r2 df
r2_df <- filtered_df %>% select(OTHER_VARIANT,R2) %>% rename(variant_id=OTHER_VARIANT,r2=R2)
r2_df <- rbind(r2_df,data.frame(variant_id=lead_variant_raw,r2=1))
r2_df <- r2_df %>% mutate(variant_id=gsub("\\:",".",variant_id))

##########################
# PREPARING FOR PLOTTING #
##########################
# assembling marginal/main/interaction plot DF
regional_plot_df<-data.frame(results)
rownames(regional_plot_df) <- NULL
colnames(regional_plot_df) <- c("BETA","SE","T","P","variant_id","phenotype_id","class")
# making an indicator variable for whether an observation is about the lead variant
regional_plot_df <- regional_plot_df %>% mutate(isleadvariant=ifelse(variant_id==lead_variant_period,1,0))
# preparing and parsing the table
regional_plot_df <- regional_plot_df %>% mutate(BETA=as.numeric(BETA),P=as.numeric(P),SE=as.numeric(SE),T=as.numeric(T),Sign = ifelse(BETA > 0, "Positive", ifelse(BETA < 0, "Negative", NA))) %>% tidyr::separate(variant_id, into = c("chr_char","pos_char","a2","a1"),sep="\\.",remove=F) %>% mutate(chr_num=as.numeric(gsub("chr","",chr_char)),pos_num=as.numeric(pos_char)) %>% mutate(
  class=factor(regional_plot_df$class, 
               levels = c("Marginal Effect", 
                          "Marginal Effect (conditional)", 
                          "Main Effect", 
                          "Main Effect (conditional)", 
                          "Interaction Effect", 
                          "Interaction Effect (conditional)"))
)

###################################
# PLOTTING THE COMPREHENSIVE PLOT #
###################################
library(ggplot2)
library(ggnewscale)
library(scales)
library(ggrepel)
library(dplyr)
library(stringr)

# Filter out conditional classes
regional_plot_df <- regional_plot_df %>%
  filter(!str_detect(tolower(class), "conditional"))

# Determine x-axis window range
x_min <- min(window_variants$POS)
x_max <- max(window_variants$POS)
x_range <- x_max - x_min

# Lead variant label data
lead_labels_df <- regional_plot_df %>%
  filter(variant_id == lead_variant_period) %>%
  mutate(
    rsid = query_rsid,
    x_point = pos_num,
    x_label = pos_num + 0.03 * x_range,
    y_coord = -log10(P)
  )

# Create the plot
plot <- ggplot() +
  
  # Plot all points with fixed sizes
  geom_point(data = regional_plot_df, 
             aes(x = pos_num / 1e6,
                 y = -log10(P),
                 color = Sign,
                 shape = factor(isleadvariant),
                 size = ifelse(isleadvariant == 1, 6, 3))) +  # ← fixed size, bigger diamond
  scale_color_manual(name = "Sign", values = c("Positive" = "blue", "Negative" = "red"), na.translate = FALSE) +
  scale_shape_manual(values = c("0" = 16, "1" = 18)) +  # Circle for 0, Diamond for 1
  scale_size_identity() +
  guides(shape = "none") +
  
  # Horizontal line to label
  geom_segment(
    data = lead_labels_df,
    aes(x = x_point / 1e6, xend = x_label / 1e6,
        y = y_coord, yend = y_coord),
    color = "black", linewidth = 0.4
  ) +
  
  # Label to the right of lead variant
  geom_label(
    data = lead_labels_df,
    aes(x = x_label / 1e6, y = y_coord, label = rsid),
    size = 3.5,
    fontface = "bold",
    fill = "white",
    color = "black",
    label.size = 0.4,
    hjust = 0
  ) +
  
  # Facet by class in one column
  facet_wrap(~ class, ncol = 1) +
  
  # Axis labels
  labs(y = "-log10(P)") +
  scale_x_continuous(
    name = paste0("Position on Chromosome ", num_chr, " (Mb)"),
    labels = comma_format(),
    breaks = c(x_min, x_max) / 1e6
  ) +
  
  # Theme
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.spacing = unit(0.5, "cm"),
    axis.text.x = element_text(hjust = c(0, 1)),
    legend.position = "none",
    strip.text = element_blank()
  ) 

# Save the plot
ggsave(
  filename = paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/regional_plots/figure2/regional_plot_", tmp_combination, "_",
                    gsub("\\:", ".", lead_variant_period), ".png"),
  plot = plot,
  width = 5,
  height = 5,
  dpi = 600
)
