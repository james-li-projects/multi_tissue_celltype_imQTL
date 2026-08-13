import sys
import pandas as pd
import os

dataset = sys.argv[1]
path_prefix = sys.argv[2]

# Set path
os.chdir(f"/gpfs/data/pierce-lab/james.li/imQTL/output/{dataset}/imQTL/full_assoc")

# Load significant pair IDs once into a set (fast lookup and low memory)
sig_pair_id_path = "/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/cpg_variant_pairs/sig_pair_id.list"
sig_pair_ids = set(pd.read_csv(sig_pair_id_path, header=None)[0])

# Prepare output file paths
sig_b_gi_path = f"/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/b_gi/sig_{dataset}_{path_prefix}.tsv"
sig_b_gi_se_path = f"/gpfs/data/pierce-lab/james.li/imQTL/output/cross_tissue_compare/b_gi_se/sig_{dataset}_{path_prefix}.tsv"

# Open output files for writing and write headers
with open(sig_b_gi_path, 'w') as b_gi_file, open(sig_b_gi_se_path, 'w') as b_gi_se_file:
    b_gi_file.write("pair_id\tb_gi\n")
    b_gi_se_file.write("pair_id\tb_gi_se\n")
    
    for i in range(1, 23):
        print(f"Processing chromosome {i}")
        parquet_path = f"{path_prefix}.chunk{i}.cis_qtl_pairs.chr{i}.parquet"
        tmp_df = pd.read_parquet(parquet_path, columns=["phenotype_id", "variant_id", "b_gi", "b_gi_se"])
        tmp_df["pair_id"] = tmp_df["phenotype_id"] + "_" + tmp_df["variant_id"]
        
        # Filter early
        tmp_df = tmp_df[tmp_df["pair_id"].isin(sig_pair_ids)]

        # Write to files in append mode
        tmp_df[["pair_id", "b_gi"]].to_csv(b_gi_file, sep="\t", header=False, index=False)
        tmp_df[["pair_id", "b_gi_se"]].to_csv(b_gi_se_file, sep="\t", header=False, index=False)
