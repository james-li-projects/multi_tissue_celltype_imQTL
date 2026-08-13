import sys
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import gc

pd.set_option('display.max_columns', None)

# Define paths
gtex_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/imQTL/top_assoc"
heals_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/imQTL/top_assoc"
full_assoc_base = "/gpfs/data/pierce-lab/james.li/imQTL/output"
plot_out_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/pleiotropy_analysis/plots"
tsv_out_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/pleiotropy_analysis/tables"

# Collect combinations
wide_parsed_imqtl = pd.read_csv("/gpfs/data/pierce-lab/james.li/imQTL/output/parsed_imqtl_effect/wide_parsed_imqtl.txt", sep="\t")
wide_parsed_imqtl["dataset_combination"] = wide_parsed_imqtl["combination"].apply(
    lambda x: "HEALS_" + x if "wb" in x else "GTEx_" + x
)
all_combinations = wide_parsed_imqtl["dataset_combination"].unique().tolist()

# Start analysis
for input_string in all_combinations:
    try:
        dataset, tissue, celltype = input_string.split("_")
        print(f"\nProcessing {dataset} - {tissue} - {celltype}")

        # Load lead variants: instead of filtering based on p-values,
        # use pre-filtered variant list from wide_parsed_imqtl
        variant_ids = wide_parsed_imqtl.loc[wide_parsed_imqtl["dataset_combination"] == input_string, "variant_id"].dropna().unique()
        imqtl_variants = set(variant_ids)
        gc.collect()
        # Collect variant–CpG pairs
        variant_cpg_pairs = []
        for i in range(1, 23):
            chr_str = str(i)
            print(f"  Processing chr{chr_str}")
            parquet_path = f"{full_assoc_base}/{dataset}/imQTL/full_assoc/tensorQTL_imQTL_{tissue}_{celltype}.chunk{chr_str}.cis_qtl_pairs.chr{chr_str}.parquet"
            try:
                chunk = pd.read_parquet(parquet_path, engine='pyarrow', columns=["phenotype_id", "variant_id", "pval_gi"])
                chunk = chunk[chunk["pval_gi"] < 1e-4]
                chunk = chunk[chunk["variant_id"].isin(imqtl_variants)]
                chunk["variant_id"] = chunk["variant_id"].astype("category")
                chunk["phenotype_id"] = chunk["phenotype_id"].astype("category")
                
                variant_cpg_pairs.extend(chunk[["variant_id", "phenotype_id"]].drop_duplicates().itertuples(index=False, name=None))
                del chunk
                gc.collect()
            except Exception as e:
                print(f"    Error reading {parquet_path}: {e}")
                continue
        # Convert to DataFrame
        variant_cpg_df = pd.DataFrame(variant_cpg_pairs, columns=["variant_id", "phenotype_id"])
        # Count CpGs per variant
        count_df = variant_cpg_df.groupby("variant_id")["phenotype_id"].nunique().reset_index()
        count_df.rename(columns={"phenotype_id": "count"}, inplace=True)
        # Plot violin plot
        plt.figure(figsize=(10, 6), dpi=300)
        sns.violinplot(x='count', data=count_df, orient='h', inner='box', color="skyblue")
        plt.title("Number of CpGs regulated by imQTL lead variants", fontsize=14)
        plt.xlabel("Count", fontsize=12)
        plt.tight_layout()
        plot_path = f"{plot_out_dir}/{dataset}_{tissue}_{celltype}_violin_plot.png"
        plt.savefig(plot_path, format='png', dpi=300)
        plt.close()
        # Save variant–phenotype pairs instead of just counts
        tsv_output_path = f"{tsv_out_dir}/{dataset}_{tissue}_{celltype}_variant_cpg_pairs.tsv"
        variant_cpg_df.to_csv(tsv_output_path, sep='\t', index=False)
        # Optionally save count table too
        count_output_path = f"{tsv_out_dir}/{dataset}_{tissue}_{celltype}_variant_counts.tsv"
        count_df.to_csv(count_output_path, sep='\t', index=False)
        del variant_cpg_df, variant_cpg_pairs, count_df
        gc.collect()
    except Exception as e:
        print(f"  ERROR processing {input_string}: {e}")
        continue
