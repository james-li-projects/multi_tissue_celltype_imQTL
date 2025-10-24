import sys
import pandas as pd
import os
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from scipy.stats import chi2

# importing arguments for dataset and path_prefix for the tissue/cell type combination 
dataset = sys.argv[1]
path_prefix = sys.argv[2]

# set path
os.chdir(f"/gpfs/data/pierce-lab/james.li/imQTL/output/{dataset}/imQTL/full_assoc")

#####################################
# Collecting all p-values efficiently
pval_list = []

for i in range(1, 23):
    print(f"Processing chromosome {i}")
    parquet_path = f"{path_prefix}.chunk{i}.cis_qtl_pairs.chr{i}.parquet"
    
    # Load only the necessary column
    tmp_df = pd.read_parquet(parquet_path, columns=["pval_gi"])
    
    # Convert to numpy array early and delete DataFrame
    pval_array = tmp_df["pval_gi"].to_numpy(dtype=np.float64, copy=False)
    del tmp_df
    
    pval_list.append(pval_array)

# Concatenate all arrays into one
all_pvalues = np.concatenate(pval_list)
del pval_list  # Free memory

#####################################
# Compute inflation factor
def calculate_inflation_factor(pvalues):
    chisq = chi2.ppf(1 - pvalues, df=1)
    median_chisq = np.median(chisq)
    theoretical_chisq = chi2.ppf(0.5, df=1)
    return median_chisq / theoretical_chisq

inflation_factor = calculate_inflation_factor(all_pvalues)
rounded_inflation_factor = f"{round(float(inflation_factor), 3):.3f}"

#####################################
# Generate Q-Q plot
# Q-Q plot setup
# Define tissue and cell type maps
tissue_map = {
    "breast": "Breast", "colon": "Colon", "lung": "Lung", "kidney": "Kidney",
    "prostate": "Prostate", "wb": "Whole Blood", "ovary": "Ovary"
}

celltype_map = {
    "Basal": "Basal epithelium", "Luminal": "Luminal epithelium", "Epi": "Epithelial cell",
    "BE": "Basal epithelium", "LE": "Luminal epithelium", "Leu": "Leukocyte", "SM": "Smooth muscle cell",
    "EC": "Endothelial cell", "Fat": "Adipocyte", "Fib": "Fibroblast", "Lym": "Lymphocyte",
    "MP": "Macrophage", "Macro": "Macrophage", "Mono": "Monocyte", "Mye": "Myeloid cell",
    "Gran": "Granulocyte", "Stromal": "Stromal cell", "Endo": "Endothelial cell",
    "B": "B cell", "CD4T": "CD4+ T Cell", "CD8T": "CD8+ T Cell", "NK": "NK cell", "Eosino": "Eosinophil",
    "Neutro": "Neutrophil", "EndoC": "Endothelial cell", "IC": "Immune cells"
}

# Q-Q plot setup
res = stats.probplot(all_pvalues, dist="uniform", plot=None, fit=False)
expected_quantiles = res[0]
observed_quantiles = res[1]

log10_expected = -np.log10(expected_quantiles)
log10_observed = -np.log10(observed_quantiles)

min_max_length = min(log10_expected.max(), log10_observed.max())
identity_line_x = np.linspace(0, min_max_length, 100)
identity_line_y = identity_line_x

plt.figure(figsize=(6, 6))
plt.scatter(log10_expected, log10_observed, s=5)
plt.plot(identity_line_x, identity_line_y, color='red', linestyle='--', label='Identity Line')

# Extract and map combination
raw_combination = parquet_path.replace("tensorQTL_imQTL_", "").replace(".chunk22.cis_qtl_pairs.chr22.parquet", "")
tissue_part, celltype_part = raw_combination.split("_", 1)
tissue_pretty = tissue_map.get(tissue_part.lower(), tissue_part)
celltype_pretty = celltype_map.get(celltype_part, celltype_part)
formatted_combination = f"{tissue_pretty} - {celltype_pretty}"

# Labels and output
plt.title(f'{formatted_combination}', fontsize=20)
plt.xlabel('Expected -log10(p-value)', fontsize=20)
plt.ylabel('Observed -log10(p-value)', fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.text(0.3, 0.8, f"λ: {rounded_inflation_factor}", 
         transform=plt.gca().transAxes, fontsize=20, color="black", ha='right')

target_path = f"/gpfs/data/pierce-lab/james.li/imQTL/output/qqplots/{dataset}_{raw_combination}.png"
plt.savefig(target_path, dpi=300)
plt.close()
