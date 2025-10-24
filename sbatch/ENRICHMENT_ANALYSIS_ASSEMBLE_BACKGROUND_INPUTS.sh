###########################################
############ BACKGROUND/TESTED ############
###########################################

############################################################
# extracting background variants tested when mapping mQTLs 
module load gcc
module load miniconda3
conda activate pytorch_env
python
import pandas as pd
import pyarrow.parquet as pq
import os
from glob import glob
# Base directories
base_input_dir = "/gpfs/data/pierce-lab/james.li/imQTL/output"
output_base = "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/intermediate"
# Dataset-specific file inclusion filters
dataset_filters = {
    "GTEx": ["colon", "lung"],
    "HEALS": ["wb"]
}
# Loop through each dataset and tissue filter
for dataset, substrings in dataset_filters.items():
    dataset_dir = os.path.join(base_input_dir, dataset, "mQTL", "full_assoc")
    all_files = glob(os.path.join(dataset_dir, "*.parquet"))
    print(f"\nProcessing dataset: {dataset} with {len(all_files)} total files")
    for substr in substrings:
        matching_files = [f for f in all_files if substr in os.path.basename(f).lower()]
        print(f"  Processing tissue: {substr} ({len(matching_files)} files)")
        variant_ids = set()
        for i, file in enumerate(matching_files, 1):
            print(f"    File {i}/{len(matching_files)}: {os.path.basename(file)}")
            try:
                df = pd.read_parquet(file, columns=["variant_id"])
                variant_ids.update(df["variant_id"].unique())
            except Exception as e:
                print(f"      Error reading {file}: {e}")
        # Save output for this dataset-tissue combo
        output_file = os.path.join(output_base, f"variant_{dataset}_{substr}.tsv")
        variant_df = pd.DataFrame(sorted(variant_ids), columns=["variant_id"])
        variant_df.to_csv(output_file, sep="\t", index=False)
        print(f"    Saved {len(variant_df)} unique variant IDs to: {output_file}")
quit()

####################################################
# creating bed files for these background variants 
conda deactivate
conda activate r_env
R
library(dplyr)
library(data.table)

# Define directory
input_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/intermediate"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed"
# Define filenames only (relative to dir_path)
filenames <- c("variant_GTEx_colon.tsv", "variant_GTEx_lung.tsv", "variant_HEALS_wb.tsv")

# Loop over files
for (fname in filenames) {
  full_path <- file.path(input_dir, fname)
  cat("\nReading file:", fname, "\n")
  dt <- fread(full_path, sep = ":", header = F)
  colnames(dt) <- c("#CHROM","POS","REF","ALT")
  variant_bed <- dt %>% mutate(ID=paste(`#CHROM`,POS,REF,ALT,sep=":"),POS2=POS) %>% select(`#CHROM`,POS,POS2,ID)
  write.table(variant_bed,file=paste0(output_dir,"/all_",gsub(".tsv",".bed",fname)),quote=F,row.names=F,col.names=F,sep="\t")
}


#####################################################
# importing and exporting bed files for tested cpgs 
tmp_dataset="GTEx"
tmp_tissue="colon"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed"
cpg_bed<-fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",tmp_dataset,"/",tmp_tissue,".bed")) %>% select(`#chr`,start,end,phenotype_id)
write.table(cpg_bed,file=paste0(output_dir,"/all_cpg_",tmp_dataset,"_",tmp_tissue,".bed"),quote=F,row.names=F,col.names=F,sep="\t")

tmp_dataset="GTEx"
tmp_tissue="lung"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed"
cpg_bed<-fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",tmp_dataset,"/",tmp_tissue,".bed")) %>% select(`#chr`,start,end,phenotype_id)
write.table(cpg_bed,file=paste0(output_dir,"/all_cpg_",tmp_dataset,"_",tmp_tissue,".bed"),quote=F,row.names=F,col.names=F,sep="\t")

tmp_dataset="HEALS"
tmp_tissue="wb"
output_dir <- "/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/ENRICHMENT_ANALYSIS/BACKGROUND/bed"
cpg_bed<-fread(paste0("/gpfs/data/pierce-lab/james.li/imQTL/input/",tmp_dataset,"/",tmp_tissue,".bed")) %>% select(`#chr`,start,end,phenotype_id)
write.table(cpg_bed,file=paste0(output_dir,"/all_cpg_",tmp_dataset,"_",tmp_tissue,".bed"),quote=F,row.names=F,col.names=F,sep="\t")
