library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpmisc)

# loading imQTL/mQTL counts
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_cpg_df.RData")
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/mQTL_feature_df.RData")
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/imQTL_cpg_df.RData")
load("/gpfs/data/pierce-lab/james.li/imQTL/output/analysis/imQTL_feature_df.RData")
imQTL_feature_df <- imQTL_feature_df %>% filter(num_imQTL!=0)

############################
all_combination_plot_df <- data.frame()
for (imQTL_combination in unique(imQTL_cpg_df$imQTL_vec)) {
  
  print(imQTL_combination)
  # identifying filenames for each imQTL and mQTL file
  imQTL_filename <- unique((imQTL_cpg_df %>% filter(imQTL_vec == imQTL_combination))$combination)
  tissue <- unique((imQTL_cpg_df %>% filter(imQTL_vec == imQTL_combination))$Tissue)
  mQTL_filename <- unique((mQTL_cpg_df %>% filter(Tissue == tissue))$file_name)
  necessary_dataset <- unique((imQTL_cpg_df %>% filter(imQTL_vec == imQTL_combination))$Dataset)
  
  # extracting mQTL associations between each CpG and lead variant
  if (necessary_dataset=="GTEx") {
    mQTL_top_assoc <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/mQTL/top_assoc/",mQTL_filename)) %>% mutate(pval_adj_bh=p.adjust(pval_beta,method="fdr")) %>% filter(pval_adj_bh<0.05) %>% mutate(midpoint = (start_distance+end_distance)/2) %>% mutate(`Distance between lead variant and CpG (kb)`=midpoint/1000) %>% mutate(QTL_type="mQTL in matching tissue",Combination=imQTL_combination)
  } else if (necessary_dataset=="HEALS") {
    mQTL_top_assoc <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/mQTL/top_assoc/",mQTL_filename)) %>% mutate(pval_adj_bh=p.adjust(pval_beta,method="fdr")) %>% filter(pval_adj_bh<0.05) %>% mutate(midpoint = (start_distance+end_distance)/2) %>% mutate(`Distance between lead variant and CpG (kb)`=midpoint/1000) %>% mutate(QTL_type="mQTL in matching tissue",Combination=imQTL_combination)
  }
  
  # extracting imQTL associations between each lead variant and CpG (lead variant - CpG)
  if (necessary_dataset=="GTEx") {
    imQTL_top_assoc <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/imQTL/top_assoc/",imQTL_filename)) %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni")) %>% filter(pval_adj_bonf<0.05) %>% mutate(midpoint = -1*(start_distance+end_distance)/2) %>% mutate(`Distance between lead variant and CpG (kb)`=midpoint/1000) %>% mutate(QTL_type="imQTL",Combination=imQTL_combination)
  } else if (necessary_dataset=="HEALS") {
    imQTL_top_assoc <- fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc/",imQTL_filename)) %>% mutate(pval_adj_bonf=p.adjust(pval_emt,method="bonferroni")) %>% filter(pval_adj_bonf<0.05) %>% mutate(midpoint = -1*(start_distance+end_distance)/2) %>% mutate(`Distance between lead variant and CpG (kb)`=midpoint/1000) %>% mutate(QTL_type="imQTL",Combination=imQTL_combination)
  }
  
  # creating histogram distances of imQTLs and mQTLs
  plot_imQTL_top_assoc <- imQTL_top_assoc %>% select(`Distance between lead variant and CpG (kb)`,QTL_type,Combination,phenotype_id)
  plot_mQTL_top_assoc <- mQTL_top_assoc %>% select(`Distance between lead variant and CpG (kb)`,QTL_type,Combination,phenotype_id)
  
  plot_density_df <- rbind(plot_imQTL_top_assoc,plot_mQTL_top_assoc)
  all_combination_plot_df <- rbind(all_combination_plot_df,plot_density_df)
  
  #png(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/DISTANCE_PLOTS/distance_density_plot_",imQTL_combination,".png"),res=600,height=3,width=6,units="in")
  #p <- ggplot(data=plot_density_df,aes(x=`Distance between lead variant and CpG (kb)`,fill=QTL_type)) + facet_wrap(~QTL_type, scales = "free_y") + xlim(-500, 500) + theme_classic() + geom_histogram() + theme(legend.title=element_blank(),legend.position = "bottom")
  #print(p)
  #dev.off()
  
  library(dplyr)
  library(ggplot2)
  
  # Parameters
  output_path <- paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/DISTANCE_PLOTS/distance_density_plot_", imQTL_combination, ".png")
  binwidth <- 20
  primary_type <- "imQTL"
  secondary_type <- "mQTL in matching tissue"
  x_limits <- c(-500, 500)
  bin_edges <- seq(x_limits[1], x_limits[2], by = binwidth)
  
  # Filter
  df_primary <- plot_density_df %>% filter(QTL_type == primary_type)
  df_secondary <- plot_density_df %>% filter(QTL_type == secondary_type)
  
  # Proceed only if both groups have enough data
  if (nrow(df_primary) >= 2 & nrow(df_secondary) >= 2) {
    
    # Create histograms manually
    hist_primary <- hist(df_primary$`Distance between lead variant and CpG (kb)`,
                         breaks = bin_edges, plot = FALSE)
    hist_secondary <- hist(df_secondary$`Distance between lead variant and CpG (kb)`,
                           breaks = bin_edges, plot = FALSE)
    
    primary_df <- data.frame(x = hist_primary$mids, y = hist_primary$counts)
    secondary_df <- data.frame(x = hist_secondary$mids, y = hist_secondary$counts)
    
    # Handle division by zero if no secondary counts
    scale_factor <- ifelse(max(secondary_df$y) > 0,
                           max(primary_df$y) / max(secondary_df$y),
                           1)
    
    # Plot
    png(output_path, res = 600, height = 3, width = 4, units = "in")
    
    p <- ggplot() +
      geom_col(data = primary_df, aes(x = x, y = y), fill = "#F8766D", alpha = 0.5, width = binwidth) +
      geom_col(data = secondary_df, aes(x = x, y = y * scale_factor), fill = "#00BFC4", alpha = 0.5, width = binwidth) +
      xlim(x_limits) +
      scale_y_continuous(
        name = "imQTL count",
        sec.axis = sec_axis(~ . / scale_factor, name = "mQTL count")
      ) +
      xlab("CpG distance to lead variant (kb)") +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.y.left = element_text(color = "#F8766D"),
        axis.title.y.right = element_text(color = "#00BFC4"),
        plot.title = element_text(hjust = 0.5)
      ) + labs(title=imQTL_combination)
    
    print(p)
    dev.off()
    
  } else {
    message("Skipping plot: Not enough data for one or both QTL types in ", imQTL_combination)
  }
}

# plotting all COMBINED imqtls vs mqtls
png(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/DISTANCE_PLOTS/distance_density_plot_COMBINED_COMBINATIONS.png"), res=1200, height=5, width=5, units="in")
p <- ggplot(data=all_combination_plot_df, aes(x=`Distance between lead variant and CpG (kb)`, fill=QTL_type)) + 
  theme_classic() + 
  geom_histogram(bins=30, position = "identity", alpha = 0.5, mapping = aes(y = stat(density))) + 
  labs(fill = "Type of QTL") +  # Set the legend title
  theme(legend.position = "bottom") + # Move legend to the bottom
  ylab("Density")
print(p)
dev.off()

# printing out summary and SD of absolute distances of imQTLs
print(summary(abs((all_combination_plot_df%>%filter(QTL_type=="imQTL"))$`Distance between lead variant and CpG (kb)`)))
print(sd(abs((all_combination_plot_df%>%filter(QTL_type=="imQTL"))$`Distance between lead variant and CpG (kb)`)))

# printing out summary and SD of absolute distances of mQTLs in matching tissues
print(summary(abs((all_combination_plot_df%>%filter(QTL_type=="mQTL in matching tissue"))$`Distance between lead variant and CpG (kb)`)))
print(sd(abs((all_combination_plot_df%>%filter(QTL_type=="mQTL in matching tissue"))$`Distance between lead variant and CpG (kb)`)))

#########################################
# creating a dataframe of all combinations
all_combination_plot_AND_ALL_COMB_df<- rbind(
  all_combination_plot_df,
  all_combination_plot_df %>% mutate(Combination="All Combinations")
) %>% mutate(QTL_type=ifelse(QTL_type=="mQTL","mQTL in matching tissue",QTL_type))

# plotting all tissue/cell type combinations
png(paste0("/gpfs/data/pierce-lab/james.li/imQTL/output/DISTANCE_PLOTS/distance_density_plot_ALL_COMBINATIONS.png"),res=600,height=8,width=9,units="in")
p <- ggplot(data=all_combination_plot_AND_ALL_COMB_df,aes(x=`Distance between lead variant and CpG (kb)`,fill=QTL_type)) + facet_wrap(~Combination, scales = "free",ncol=4) + theme_classic() + geom_histogram(bins=30, position = "identity", alpha = 0.5, mapping = aes(y = stat(density))) + theme(legend.title=element_blank()) + ylab("Density") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # Add borders
  )
print(p)
dev.off()

###########################
# mean absolute distance barplots
library(ggplot2)
library(dplyr)
summary_df <- all_combination_plot_AND_ALL_COMB_df %>% rename(`CpG distance from lead variant (kb)`=`Distance between lead variant and CpG (kb)`) %>%
  mutate(abs_distance = abs(`CpG distance from lead variant (kb)`)) %>%
  group_by(Combination, QTL_type) %>%
  summarize(
    mean_dist = mean(abs_distance, na.rm = TRUE),
    .groups = "drop"
  )
png("/gpfs/data/pierce-lab/james.li/imQTL/output/DISTANCE_PLOTS/mean_abs_distance_barplot.png", res = 600, height = 4.25, width = 10, units = "in")
p <- ggplot(summary_df, aes(x = Combination, y = mean_dist, fill = QTL_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, alpha = 0.5) +
  geom_text(aes(label = round(mean_dist, 1)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 1.75) +
  theme_classic() +
  ylab("CpG distance to lead variant (kb)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.title = element_blank(),
        legend.position = "bottom",

strip.background = element_blank(),
strip.text = element_text(face = "bold", size = 8),
panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)
print(p)
dev.off()






#########################################
#### EXTRA CODE THAT'S NOT NECESSARY ####
#########################################
# identifying # CpGs that uniquely mapped to imQTLs and not mQTLs for each tissue/celltype combination
unique_count_df <- data.frame()
for (tmp_Combination in unique(all_combination_plot_df$Combination)) {
  print(tmp_Combination)
  tmp_Combination_df <- all_combination_plot_df %>% filter(Combination==tmp_Combination)
  cpg_list_mqtl <- unique((tmp_Combination_df %>% filter(QTL_type=="mQTL")%>%select(phenotype_id)%>%unique())$phenotype_id)
  cpg_list_imqtl <- unique((tmp_Combination_df %>% filter(QTL_type=="imQTL")%>%select(phenotype_id)%>%unique())$phenotype_id)
  only_imqtl_cpg_list <- setdiff(cpg_list_imqtl,cpg_list_mqtl)
  tmp_unique_count_df <- data.frame(
    t_ct_combination = tmp_Combination,
    unique_imqtl=length(only_imqtl_cpg_list),
    all_imqtl=length(cpg_list_imqtl),
    perc_unique_imqtl=length(only_imqtl_cpg_list)/length(cpg_list_imqtl)
  )
  unique_count_df <- rbind(unique_count_df,tmp_unique_count_df)
}


# identifying # CpGs that uniquely mapped to imQTLs and not mQTLs for each tissue
modified_all_combination_plot_df <- all_combination_plot_df %>% mutate(Combination=gsub(" - ","-",Combination)) %>% separate(Combination,into=c("Tissue","Cell Type"),sep="\\-",remove=F)
unique_count_df <- data.frame()
for (tmp_Tissue in unique(modified_all_combination_plot_df$Tissue)) {
  print(tmp_Tissue)
  tmp_Tissue_df <- modified_all_combination_plot_df %>% filter(Tissue==tmp_Tissue)
  cpg_list_mqtl <- unique((tmp_Tissue_df %>% filter(QTL_type=="mQTL")%>%select(phenotype_id)%>%unique())$phenotype_id)
  cpg_list_imqtl <- unique((tmp_Tissue_df %>% filter(QTL_type=="imQTL")%>%select(phenotype_id)%>%unique())$phenotype_id)
  only_imqtl_cpg_list <- setdiff(cpg_list_imqtl,cpg_list_mqtl)
  tmp_unique_count_df <- data.frame(
    Tissue = tmp_Tissue,
    unique_imqtl=length(only_imqtl_cpg_list),
    all_imqtl=length(cpg_list_imqtl),
    perc_unique_imqtl=length(only_imqtl_cpg_list)/length(cpg_list_imqtl)
  )
  unique_count_df <- rbind(unique_count_df,tmp_unique_count_df)
}

