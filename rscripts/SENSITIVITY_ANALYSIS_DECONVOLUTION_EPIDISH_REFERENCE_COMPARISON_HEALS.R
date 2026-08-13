library(data.table)
library(tidyverse)
library(EpiSCORE)
library(reshape)
library(EpiDISH)

meth_input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/HEALS/methylation"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac"

# loading function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
tissue_list <- c("wb")

tissue_mref_list <- vector(length=length(tissue_list),mode="list")
load("/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/EPIDISH_12CT_EPIC_PANEL/cent12CT.m.rda")
tissue_mref_list[[1]]<-cent12CT.m

# estimating cell type fractions for whole blood
for (i in 1) {
  current_tissue <- tissue_list[i]
  current_tissue_mref <- tissue_mref_list[[i]]
  
  # loading sample list for a given tissue
  mQTL_sample_list <- loadRData(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/HEALS/mQTL_sample_list_",current_tissue,".RData"))
  print(paste("Inferring cell type fractions for",current_tissue))
  beta_mat <- loadRData(paste0(meth_input_dir,"/combined_beta_modid.RData"))
  beta_mat <- beta_mat[,colnames(beta_mat)%in%mQTL_sample_list]
  wRPC.o <- epidish(beta.m = beta_mat, ref.m = current_tissue_mref, method = 'RPC')
  estF.o_cent12CT <- wRPC.o$estF
  
}

# loading previous estimates
load("/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac/estF.wb.RData")

# converting to dataframes
estF.o_cent12CT <- data.frame(estF.o_cent12CT)
estF.o <- data.frame(estF.o)

# Create collapsed dataframe
estF.o_cent12CT_collapsed <- estF.o_cent12CT %>%
  dplyr::transmute(
    B      = Bnv + Bmem,
    NK     = NK,
    CD4T   = CD4Tnv + CD4Tmem + Treg,
    CD8T   = CD8Tnv + CD8Tmem,
    Mono   = Mono,
    Neutro = Neu,
    Eosino = Eos
  )

# Preserve rownames
rownames(estF.o_cent12CT_collapsed) <- rownames(estF.o_cent12CT)

# Reorder columns to exactly match estF.o
estF.o_cent12CT_collapsed <- estF.o_cent12CT_collapsed[, colnames(estF.o)]

# Check
head(estF.o_cent12CT_collapsed, 1)





# Ensure same column order
common_ct <- intersect(colnames(estF.o), colnames(estF.o_cent12CT_collapsed))

# Ensure same row order
common_samples <- intersect(rownames(estF.o), rownames(estF.o_cent12CT_collapsed))

est1 <- estF.o[common_samples, common_ct]
est2 <- estF.o_cent12CT_collapsed[common_samples, common_ct]

# Compute correlations per cell type
cor_results <- lapply(common_ct, function(ct) {
  
  test <- cor.test(est1[, ct], est2[, ct], method = "pearson")
  
  data.frame(
    cell_type = ct,
    r  = unname(test$estimate),
    r2 = unname(test$estimate)^2,
    p_value = test$p.value
  )
})

cor_results <- do.call(rbind, cor_results)

cor_results

############################
# plotting code
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---- align rows/cols ----
common_ct <- intersect(colnames(estF.o), colnames(estF.o_cent12CT_collapsed))
common_samples <- intersect(rownames(estF.o), rownames(estF.o_cent12CT_collapsed))

est1 <- estF.o[common_samples, common_ct, drop = FALSE]
est2 <- estF.o_cent12CT_collapsed[common_samples, common_ct, drop = FALSE]

# ---- long format ----
df <- tibble(
  sample = common_samples,
  !!!as.data.frame(est1)
) %>%
  pivot_longer(cols = all_of(common_ct), names_to = "cell_type", values_to = "x") %>%
  left_join(
    tibble(
      sample = common_samples,
      !!!as.data.frame(est2)
    ) %>% pivot_longer(cols = all_of(common_ct), names_to = "cell_type", values_to = "y"),
    by = c("sample", "cell_type")
  )

# ---- FORCE Neutro first ----
ct_order <- c("Neutro", setdiff(common_ct, "Neutro"))
df$cell_type <- factor(df$cell_type, levels = ct_order)

# ---- stats per facet ----
stats <- df %>%
  group_by(cell_type) %>%
  summarise(
    r = cor(x, y, use = "pairwise.complete.obs"),
    p = {
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) >= 3) cor.test(x[ok], y[ok])$p.value else NA_real_
    },
    x_min = min(x, na.rm = TRUE),
    x_max = max(x, na.rm = TRUE),
    y_min = min(y, na.rm = TRUE),
    y_max = max(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "r=", sprintf("%.3f", r),
      "\n",
      "p=", ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))
    ),
    x_lab = x_min + 0.03 * (x_max - x_min),
    y_lab = y_max - 0.03 * (y_max - y_min)
  )

# ---- plot ----
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  geom_text(
    data = stats,
    aes(x = x_lab, y = y_lab, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 4
  ) +
  facet_wrap(~ cell_type, scales = "free", nrow = 1) +
  labs(
    x = "Estimated Cell Proportion (7-cell DHS reference)",
    y = "Estimated Cell Proportion\n(12-cell immune reference,\ncollapsed to 7-cell types)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# ---- save ----
ggsave(
  filename = "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/cell_type_frac/Reinius_vs_cent12CT.png",
  plot = p,
  width = 28,
  height = 4.5,
  units = "in",
  dpi = 400,
  bg = "white"
)
