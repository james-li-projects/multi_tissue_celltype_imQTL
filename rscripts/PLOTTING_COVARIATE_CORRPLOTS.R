############################################
################### GTEX ###################
############################################

library(data.table)
library(dplyr)
library(stringr)
library(RNOmni)

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx"

tissue_list <- c("kidney","breast","colon","lung","prostate","wb","ovary")
for (current_tissue in tissue_list) {
  print(paste("PROCESSING COVARIATE FILE FOR TISSUE:",current_tissue))
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/mQTL_sample_list_",current_tissue,".RData"))
  
  # importing GTEx covariates
  current_tissue_cov <- fread("/gpfs/data/pierce-lab/GTEx/GTEx_data/Metadata/phs000424.v10.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")
  colnames(current_tissue_cov) <- as.character(current_tissue_cov[1, ])
  current_tissue_cov <- current_tissue_cov[-1, ]
  current_tissue_cov <- current_tissue_cov %>% rename(CollaboratorParticipantID=SUBJID) %>% select(CollaboratorParticipantID,SEX,AGE,BMI,MHSMKSTS,RACE) %>% filter(CollaboratorParticipantID%in%mQTL_sample_list) %>% mutate(SEX=as.numeric(SEX)) %>% mutate(SEX = SEX-1)
  
  # importing genotyping PCs
  genotyping_PC <- fread(paste0(input_dir,"/genetic_data/PLINK2_PCA_30.eigenvec"),header=T) %>% select(IID, PC1, PC2, PC3, PC4, PC5)
  colnames(genotyping_PC)[1:6] <- c("CollaboratorParticipantID","G_PC1","G_PC2","G_PC3","G_PC4","G_PC5")
  genotyping_PC <- genotyping_PC %>% select(CollaboratorParticipantID,G_PC1,G_PC2,G_PC3,G_PC4,G_PC5)
  # importing tissue-specific PCs
  current_tissue_PC <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/QTL_PCS_",current_tissue,".RData"))
  current_tissue_PC_df <- data.frame(current_tissue_PC)
  current_tissue_PC_df$CollaboratorParticipantID <- rownames(current_tissue_PC_df)
  # assembling and outputing covariate files
  combined_cov <- inner_join(inner_join(current_tissue_cov,genotyping_PC,by=c("CollaboratorParticipantID")),current_tissue_PC_df,by=c("CollaboratorParticipantID"))
  
  # keeping the order of the intersecting sample list
  combined_cov <- data.frame(combined_cov)
  rownames(combined_cov) <- combined_cov$CollaboratorParticipantID
  combined_cov <- combined_cov[mQTL_sample_list,]
  
  # importing and joining cell types
  load(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac/estF.",current_tissue,".RData"))
  estF.o <- data.frame(estF.o)
  estF.o$CollaboratorParticipantID <- rownames(estF.o)
  combined_cov <- combined_cov %>% inner_join(estF.o,by=c("CollaboratorParticipantID"))
  
  # updating column names of celltypes, G_PCs, and PCs
  celltype_map <- c(
    "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
    "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
    "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
    "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
    "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
    "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
    "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
  )
  
  # Original column names
  colnames_new <- colnames(combined_cov)
  
  # Step 1: Rename G_PC# to Genotyping PC#
  colnames_new <- gsub("^G_PC(\\d+)$", "Genotyping PC\\1", colnames_new)
  
  # Step 2: Rename PC# to DNAm PC#, but avoid already-changed Genotyping PCs
  colnames_new <- gsub("(?<!Genotyping )PC(\\d+)$", "DNAm PC\\1", colnames_new, perl = TRUE)
  
  # Step 3: Find last DNAm PC column to determine where cell type columns begin
  last_dnampc_index <- max(grep("^DNAm PC\\d+$", colnames_new))
  
  # Step 4: Rename remaining columns based on celltype_map
  for (i in (last_dnampc_index + 1):length(colnames_new)) {
    if (colnames_new[i] %in% names(celltype_map)) {
      colnames_new[i] <- celltype_map[[colnames_new[i]]]
    }
  }
  
  # Assign updated column names
  colnames(combined_cov) <- colnames_new
  
  # modifying smoking variable
  combined_cov <- combined_cov %>% mutate(MHSMKSTS=ifelse(MHSMKSTS=="Yes",1,ifelse(MHSMKSTS=="No",0,NA))) %>% rename(`Smoking status`=MHSMKSTS) 
  
  # if the tissue is prostate or ovary or breast, exclude the sex variable
  if (current_tissue %in% c("prostate","ovary","breast")) {
    combined_cov <- combined_cov %>% select(-SEX)
  } else {
    print("Tissue does not only have samples that are sex-specific (e.g. prostate, ovary, or breast)")
    # modifying Sex column names
    combined_cov <- combined_cov %>% rename(Sex=SEX)
    
  }
  
  # modifying other covariate column names
  combined_cov <- combined_cov %>% rename(Age=AGE,Race=RACE)
  
  # making correlation plot
  library(corrplot)
  
  # Remove ID column and genotyping PCs
  cov_numeric <- combined_cov[, !grepl("^Genotyping|^CollaboratorParticipantID$", colnames(combined_cov))]
  
  # printing out summary of particpant characteristics
  numeric_converted <- data.frame(lapply(
    cov_numeric[, 1:(which(names(cov_numeric) == "Race") - 1)],
    function(x) suppressWarnings(as.numeric(as.character(x)))
  ))
  summary_with_sd <- t(sapply(numeric_converted, function(x) {
    stats <- summary(x)[1:6]  # Min., 1st Qu., Median, Mean, 3rd Qu., Max.
    sd_val <- sd(x, na.rm = TRUE)
    c(stats, SD = sd_val)
  }))
  summary_df <- as.data.frame(summary_with_sd)
  summary_df <- round(summary_df, 2)
  print(summary_df)
  cat("\nRace distribution:\n")
  print(table(cov_numeric$Race, useNA = "ifany"))
  if ("Sex" %in% colnames(cov_numeric)) {
    cat("\nSex distribution:\n")
    print(table(cov_numeric$Sex, useNA = "ifany"))
  }
  cat("\nSmoking distribution:\n")
  print(table(cov_numeric$`Smoking status`, useNA = "ifany"))
  
  
  # Ensure numeric conversion (coerce categorical like "Yes"/"No" to NA to avoid errors)
  cov_numeric[] <- lapply(cov_numeric, function(x) {
    if (is.character(x) || is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else x
  })
  
  # Identify DNAm PC columns
  dnam_pc_cols <- grep("^DNAm PC", colnames(cov_numeric), value = TRUE)
  
  # Identify all other columns (rows in the plot)
  other_cols <- setdiff(colnames(cov_numeric), dnam_pc_cols)
  
  # Subset data
  dnam_pc_matrix <- as.matrix(cov_numeric[, dnam_pc_cols])
  other_matrix <- as.matrix(cov_numeric[, other_cols])
  
  # Compute correlation
  cor_matrix <- cor(other_matrix, dnam_pc_matrix, use = "pairwise.complete.obs")
  
  # Output file path
  out_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/covariate_corrplot/GTEX_DNAmPC_covariate_corrplot_", current_tissue, ".png")
  
  # Save corrplot as high-quality PNG
  png(out_path, width = 2000, height = 1800, res = 300)
  corrplot(
    cor_matrix,
    is.corr = TRUE,
    method = "color", 
    tl.col = "black",
    tl.cex = 0.8,
    tl.srt = 45, 
    cl.cex = 0.8, 
    addCoef.col = "black",   
    number.cex = 0.7     
  )
  dev.off()
}


##############################################
################### HEALS ###################
##############################################
library(data.table)
library(dplyr)
library(stringr)
library(RNOmni)
library(tidyr)

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS"

tissue_list <- c("wb")
for (current_tissue in tissue_list) {
  print(paste("PROCESSING COVARIATE FILE FOR TISSUE:",current_tissue))
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/mQTL_sample_list_",current_tissue,".RData"))
  
  # importing HEALS covariates
  current_tissue_cov <- loadRData(paste0(input_dir,"/covariates/combined_cov_modid.RData"))
  current_tissue_cov <- current_tissue_cov %>% filter(Sample_Name%in%mQTL_sample_list) %>% rename(CollaboratorParticipantID=Sample_Name,AGE=Age,SEX=Sex) %>% separate(CollaboratorParticipantID,into=c("id1","Sample_Name"),sep="_",remove=F) %>% mutate(Sample_Name=paste0("X",Sample_Name))
  
  # importing BMI
  load("/gpfs/data/phs/groups/Projects/GEMS/Zheng/data/combined_cov_BMI.RData")
  combined_cov_BMI <- combined_cov_BMI %>% select(Sample_Name,bmi) %>% rename(BMI=bmi)
  
  # joining BMI to other covariates columns
  current_tissue_cov <- current_tissue_cov %>% inner_join(combined_cov_BMI,by=c("Sample_Name"))
  
  # changing smoking status to a binary variable
  current_tissue_cov <- current_tissue_cov %>% mutate(MHSMKSTS=ifelse(cigsmoke %in% c(1,2),1,ifelse(cigsmoke == 0, 0, NA)))
  
  # making a dummy race variable
  current_tissue_cov <- current_tissue_cov %>% mutate(RACE="SAS")
  
  # selecting appropriate columns
  current_tissue_cov <- current_tissue_cov %>% select(CollaboratorParticipantID,SEX,AGE,BMI,MHSMKSTS,RACE)
  
  # importing genotyping PCs
  genotyping_PC <- fread(paste0(input_dir,"/genetic_data/PLINK2_PCA_30.eigenvec"),header=T) %>% select(IID, PC1, PC2, PC3, PC4, PC5)
  colnames(genotyping_PC)[1:6] <- c("CollaboratorParticipantID","G_PC1","G_PC2","G_PC3","G_PC4","G_PC5")
  genotyping_PC <- genotyping_PC %>% select(CollaboratorParticipantID,G_PC1,G_PC2,G_PC3,G_PC4,G_PC5)
  # importing tissue-specific PCs
  current_tissue_PC <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/QTL_PCS_",current_tissue,".RData"))
  current_tissue_PC_df <- data.frame(current_tissue_PC)
  current_tissue_PC_df$CollaboratorParticipantID <- rownames(current_tissue_PC_df)
  # assembling and outputing covariate files
  combined_cov <- inner_join(inner_join(current_tissue_cov,genotyping_PC,by=c("CollaboratorParticipantID")),current_tissue_PC_df,by=c("CollaboratorParticipantID"))
  
  # keeping the order of the intersecting sample list
  combined_cov <- data.frame(combined_cov)
  rownames(combined_cov) <- combined_cov$CollaboratorParticipantID
  combined_cov <- combined_cov[mQTL_sample_list,]
  
  # importing and joining cell types
  load(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac/estF.",current_tissue,".RData"))
  estF.o <- data.frame(estF.o)
  estF.o$CollaboratorParticipantID <- rownames(estF.o)
  combined_cov <- combined_cov %>% inner_join(estF.o,by=c("CollaboratorParticipantID"))
  
  # updating column names of celltypes, G_PCs, and PCs
  celltype_map <- c(
    "Basal" = "Basal epithelium", "Luminal" = "Luminal epithelium", "Epi" = "Epithelial cell",
    "BE" = "Basal epithelium", "LE" = "Luminal epithelium", "Leu" = "Leukocyte", "SM" = "Smooth muscle cell",
    "EC" = "Endothelial cell", "Fat" = "Adipocyte", "Fib" = "Fibroblast", "Lym" = "Lymphocyte",
    "MP" = "Macrophage", "Macro" = "Macrophage", "Mono" = "Monocyte", "Mye" = "Myeloid cell",
    "Gran" = "Granulocyte", "Stromal" = "Stromal cell", "Endo" = "Endothelial cell",
    "B" = "B cell", "CD4T" = "CD4+ T Cell", "CD8T" = "CD8+ T Cell", "NK" = "NK cell", "Eosino" = "Eosinophil",
    "Neutro" = "Neutrophil", "EndoC" = "Endothelial cell", "IC" = "Immune cells"
  )
  
  # Original column names
  colnames_new <- colnames(combined_cov)
  
  # Step 1: Rename G_PC# to Genotyping PC#
  colnames_new <- gsub("^G_PC(\\d+)$", "Genotyping PC\\1", colnames_new)
  
  # Step 2: Rename PC# to DNAm PC#, but avoid already-changed Genotyping PCs
  colnames_new <- gsub("(?<!Genotyping )PC(\\d+)$", "DNAm PC\\1", colnames_new, perl = TRUE)
  
  # Step 3: Find last DNAm PC column to determine where cell type columns begin
  last_dnampc_index <- max(grep("^DNAm PC\\d+$", colnames_new))
  
  # Step 4: Rename remaining columns based on celltype_map
  for (i in (last_dnampc_index + 1):length(colnames_new)) {
    if (colnames_new[i] %in% names(celltype_map)) {
      colnames_new[i] <- celltype_map[[colnames_new[i]]]
    }
  }
  
  # Assign updated column names
  colnames(combined_cov) <- colnames_new
  
  # modifying smoking variable
  combined_cov <- combined_cov %>% rename(`Smoking status`=MHSMKSTS) 
  
  # if the tissue is prostate or ovary or breast, exclude the sex variable
  if (current_tissue %in% c("prostate","ovary","breast")) {
    combined_cov <- combined_cov %>% select(-SEX)
  } else {
    print("Tissue does not only have samples that are sex-specific (e.g. prostate, ovary, or breast)")
    # modifying Sex column names
    combined_cov <- combined_cov %>% rename(Sex=SEX)
    
  }
  
  # modifying other covariate column names
  combined_cov <- combined_cov %>% rename(Age=AGE,Race=RACE)
  
  # making correlation plot
  library(corrplot)
  
  # Remove ID column, genotyping PCs, and Race
  cov_numeric <- combined_cov[, !grepl("^Genotyping|^CollaboratorParticipantID$|^Race$", colnames(combined_cov))]
  
  # printing out summary of particpant characteristics
  custom_summary <- function(x) {
    stats <- summary(x)
    sd_val <- sd(x, na.rm = TRUE)
    c(stats, SD = sd_val)
  }
  numeric_converted <- data.frame(lapply(cov_numeric[, 1:(which(names(cov_numeric) == "Smoking status") - 1)],
                                         function(x) suppressWarnings(as.numeric(as.character(x)))))
  summary_with_sd <- sapply(numeric_converted, custom_summary)
  print(summary_with_sd)
  
  # Ensure numeric conversion (coerce categorical like "Yes"/"No" to NA to avoid errors)
  cov_numeric[] <- lapply(cov_numeric, function(x) {
    if (is.character(x) || is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else x
  })
  cat("\nRace distribution:\n")
  print(table(cov_numeric$Race, useNA = "ifany"))
  if ("Sex" %in% colnames(cov_numeric)) {
    cat("\nSex distribution:\n")
    print(table(cov_numeric$Sex, useNA = "ifany"))
  }
  cat("\nSmoking distribution:\n")
  print(table(cov_numeric$`Smoking status`, useNA = "ifany"))
  
  # Identify DNAm PC columns
  dnam_pc_cols <- grep("^DNAm PC", colnames(cov_numeric), value = TRUE)
  
  # Identify all other columns (rows in the plot)
  other_cols <- setdiff(colnames(cov_numeric), dnam_pc_cols)
  
  # Subset data
  dnam_pc_matrix <- as.matrix(cov_numeric[, dnam_pc_cols])
  other_matrix <- as.matrix(cov_numeric[, other_cols])
  
  # Compute correlation
  cor_matrix <- cor(other_matrix, dnam_pc_matrix, use = "pairwise.complete.obs")
  
  # Output file path
  out_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/covariate_corrplot/GTEX_DNAmPC_covariate_corrplot_", current_tissue, ".png")
  
  # Save corrplot as high-quality PNG
  png(out_path, width = 2000, height = 1800, res = 300)
  corrplot(
    cor_matrix,
    is.corr = TRUE,
    method = "color", 
    tl.col = "black",
    tl.cex = 0.8,
    tl.srt = 45, 
    cl.cex = 0.8, 
    addCoef.col = "black",   
    number.cex = 0.7     
  )
  dev.off()
}

