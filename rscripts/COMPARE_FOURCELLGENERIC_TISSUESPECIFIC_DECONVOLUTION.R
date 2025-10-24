library(data.table)
library(tidyverse)
library(EpiSCORE)
library(reshape)
library(EpiDISH)

meth_input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/methylation"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac"

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

tissue_list <- c("colon","lung")

# loading centEpiFibEndoIC.m
load("/gpfs/data/pierce-lab/james.li/imQTL/teschendorff_reference/centEpiFibEndoIC.Rd")
tissue_mref_list <- vector(length=length(tissue_list)+1,mode="list")
tissue_mref_list[[1]]<-EpiSCORE::Colon_Mref.m
tissue_mref_list[[2]]<-EpiSCORE::mrefLung.m
tissue_mref_list[[3]]<-centEpiFibEndoIC.m

# estimating cell type fractions for all tissues using EpiSCORE except whole blood and ovary
for (i in 1:2) {
  # identifying current tissue type
  current_tissue <- tissue_list[i]
  current_tissue_mref <- tissue_mref_list[[i]]
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/GTEx/mQTL_sample_list_",current_tissue,".RData"))
  
  print(paste("Inferring cell type fractions for",current_tissue))
  # restrict samples 
  noob_final_BMIQ <- loadRData(paste0(meth_input_dir,"/noob_final_BMIQ_",current_tissue,"_2-6-2021.RData"))
  colnames(noob_final_BMIQ)=str_extract(colnames(noob_final_BMIQ),'GTEX-\\w+')
  noob_final_BMIQ <- noob_final_BMIQ[,colnames(noob_final_BMIQ)%in%mQTL_sample_list]
  
  
  #########################
  # TISSUE SPECIFIC PANEL #
  #########################
  # inferring cell type fractions
  constAvBetaTSS.o <-constAvBetaTSS(noob_final_BMIQ,type='850k')
  wRPC.o <- wRPC(constAvBetaTSS.o,ref=current_tissue_mref,useW=TRUE,wth=0.4,maxit=200)
  estF.o <- wRPC.o$estF
  #melt data for plotting
  difficulty_data <- estF.o %>%
    t() %>%
    as.data.frame %>% 
    t() 
  plot_df_1 <- as.data.frame.table(difficulty_data)
  colnames(plot_df) <- c("Sample_ID","Cell Type","Proportion")
  
  
  #################
  # GENERIC PANEL #
  #################
  current_tissue_mref_2 <- tissue_mref_list[[length(tissue_list)+1]]
  wRPC.o_2 <- epidish(beta.m = noob_final_BMIQ, ref.m = current_tissue_mref_2, method = 'RPC')
  estF.o_2 <- wRPC.o_2$estF
  #melt data for plotting
  difficulty_data_2 <- estF.o_2 %>%
    t() %>%
    as.data.frame %>% 
    t() 
  plot_df_2 <- as.data.frame.table(difficulty_data_2)
  
  
  #################################
  # PLOTTING AND CORRELATION CODE #
  #################################
  library(dplyr)
  library(ggplot2)
  library(tools)  # for toTitleCase()
  
  # Pretty cell type labels
  celltype_labels <- c(
    EC = "Endothelial cells",
    Epi = "Epithelial cells",
    IC = "Immune cells",
    Stromal = "Stromal cells/Fibroblasts"
  )
  
  # Clean and harmonize data
  if (i == 1) {
    plot_df_1_clean <- plot_df_1 %>%
      mutate(Var2 = as.character(Var2),
             Var2 = ifelse(Var2 %in% c("Lym", "Mye"), "IC", Var2)) %>%
      group_by(Var1, Var2) %>%
      summarize(Freq = sum(Freq), .groups = "drop")
    
    plot_df_2_clean <- plot_df_2 %>%
      mutate(Var2 = as.character(Var2),
             Var2 = recode(Var2, "EndoC" = "EC", "Fib" = "Stromal"))
    
  } else if (i == 2) {
    plot_df_1_clean <- plot_df_1 %>%
      mutate(Var2 = as.character(Var2),
             Var2 = recode(Var2, "Endo" = "EC"),
             Var2 = ifelse(Var2 %in% c("Gran", "Lym", "Macro", "Mono"), "IC", Var2)) %>%
      group_by(Var1, Var2) %>%
      summarize(Freq = sum(Freq), .groups = "drop")
    
    plot_df_2_clean <- plot_df_2 %>%
      mutate(Var2 = as.character(Var2),
             Var2 = recode(Var2, "EndoC" = "EC", "Fib" = "Stromal"))
  }
  
  # Merge
  merged_df <- inner_join(plot_df_1_clean, plot_df_2_clean,
                          by = c("Var1", "Var2"),
                          suffix = c("_1", "_2"))
  
  # Loop over each cell type
  for (ct in unique(merged_df$Var2)) {
    df_ct <- merged_df %>% filter(Var2 == ct)
    lm_model <- lm(Freq_1 ~ Freq_2, data = df_ct)
    r2_val <- round(summary(lm_model)$r.squared, 3)
    p_val <- signif(summary(lm_model)$coefficients[2, 4], 3)
    cell_label <- celltype_labels[[ct]]
    
    # Generate plot
    p <- ggplot(df_ct, aes(x = Freq_2, y = Freq_1)) +
      geom_point(size = 1.5) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      theme_classic() +
      labs(
        title = paste0(toTitleCase(tolower(current_tissue)), " - ", cell_label),
        x = paste0(cell_label, " estimated with generic four cell-type panel"),
        y = paste0(cell_label, " estimated with tissue-specific panel")
      ) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 2,
               label = paste0("r^2 == ", r2_val, " ~~~ p == ", p_val),
               parse = TRUE, size = 4)
    
    # Save
    out_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/cell_type_frac/fourcell_vs_tissuespecific/scatterplot_",
                       current_tissue, "_", ct, ".png")
    ggsave(out_path, p, width = 6, height = 6, dpi = 300)
  }
}
