#######################################
import os
import pandas as pd
from glob import glob

# Define input and output directories
input_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/mQTL/full_assoc/"
output_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/GTEx/mQTL/full_assoc_tsv/"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Define tissue keywords to filter
tissue_keywords = ["colon", "lung", "ovary", "prostate"]

# Find all .parquet files and filter for those containing the tissue keywords
parquet_files = [
    f for f in glob(os.path.join(input_dir, "*.parquet"))
    if any(tissue in os.path.basename(f).lower() for tissue in tissue_keywords)
]

# Loop through each parquet file and convert to TSV
for parquet_file in parquet_files:
    try:
        df = pd.read_parquet(parquet_file)
        base_filename = os.path.basename(parquet_file).replace(".parquet", ".tsv")
        output_path = os.path.join(output_dir, base_filename)
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Converted: {parquet_file} -> {output_path}")
    except Exception as e:
        print(f"Failed to process {parquet_file}: {e}")


#######################################
import os
import pandas as pd
from glob import glob

# Define input and output directories
input_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/mQTL/full_assoc/"
output_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/HEALS/mQTL/full_assoc_tsv/"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Find all .parquet files
parquet_files = glob(os.path.join(input_dir, "*.parquet"))

# Loop through each parquet file and convert to TSV
for parquet_file in parquet_files:
    try:
        df = pd.read_parquet(parquet_file)
        base_filename = os.path.basename(parquet_file).replace(".parquet", ".tsv")
        output_path = os.path.join(output_dir, base_filename)
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Converted: {parquet_file} -> {output_path}")
    except Exception as e:
        print(f"Failed to process {parquet_file}: {e}")
