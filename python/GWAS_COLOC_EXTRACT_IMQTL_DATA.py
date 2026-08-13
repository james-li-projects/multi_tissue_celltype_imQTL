import pandas as pd
import os
import argparse

# Set up command-line argument parser
parser = argparse.ArgumentParser(description="Extract imQTL region around phenotype and variant.")
parser.add_argument("tsv_path", type=str, help="Full path to the input TSV file")
args = parser.parse_args()

# Load the TSV file directly from the full path
df = pd.read_csv(args.tsv_path, sep='\t')

# Loop over each row
for i in range(len(df)):
    row = df.iloc[i]
    # Extract relevant values
    parquet_path = row['input_parquet_filename']
    phenotype_id = row['phenotype_id']
    pos_target = int(row['pos'])
    parsed_combination = row['parsed_combination']
    variant_id = row['variant_id']
    try:
        # Read the parquet file
        parquet_df = pd.read_parquet(parquet_path)
        # Filter for matching phenotype_id
        filtered_df = parquet_df[parquet_df['phenotype_id'] == phenotype_id].copy()
        # Split variant_id into CHR, POS, REF, ALT
        variant_parts = filtered_df['variant_id'].str.split(':', expand=True)
        filtered_df['CHR'] = variant_parts[0]
        filtered_df['POS'] = pd.to_numeric(variant_parts[1])
        filtered_df['REF'] = variant_parts[2]
        filtered_df['ALT'] = variant_parts[3]
        # Filter for POS within ±250,000 of target pos
        final_extracted_df = filtered_df[
            (filtered_df['POS'] >= pos_target - 250000) &
            (filtered_df['POS'] <= pos_target + 250000)
        ]
        # Define output path and filename
        safe_variant_id = variant_id.replace(":", "_")  # Sanitize for file name
        output_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output/GWAS_coloc/extract_imqtl_data/"
        output_path = os.path.join(output_dir, f"{parsed_combination}_{phenotype_id}_{safe_variant_id}.tsv")
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        # Save the result
        final_extracted_df.to_csv(output_path, sep='\t', index=False)
        print(f"Saved: {output_path}")
    except Exception as e:
        print(f"Error processing row {i}: {e}")

