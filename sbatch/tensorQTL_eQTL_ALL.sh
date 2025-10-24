#!/bin/bash
#SBATCH --job-name=tensorQTL_eQTL_ALL
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/tensorQTL_eQTL_ALL_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/tensorQTL_eQTL_ALL_%A.err
#SBATCH --time=24:00:00
#SBATCH -p gpuq
#SBATCH --gres gpu:1
#SBATCH --mem-per-gpu=100gb

# loading modules 
module load gcc/12.1.0
module load python/3.10.5

# importing tissue and cell type arguments
tissue=${ARGS1}

# genetic data directory
genetic_data_dir=/gpfs/data/pierce-lab/james.li/imQTL/data/GTEx/genetic_data

# covariate/methylation input directory
input_dir=/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GTEx_eQTL_v8/input

# identifying input files and paths
plink_file=${genetic_data_dir}/processed_genetic_data_chrprefix
expr_bed=${input_dir}/${tissue}.v8.normalized_expression.bed.gz
covariates_file=${input_dir}/${tissue}.v8.covariates.txt

# setting output file prefix
prefix="/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/GTEx_eQTL_v8/output/parquet/tensorQTL_eQTL_${tissue}"

# mapping imQTLs
python3 -m tensorqtl ${plink_file} ${expr_bed} ${prefix} \
    --covariates ${covariates_file} \
    --window 1000000 \
    --mode cis_nominal \
    --maf_threshold 0.01 \
    --seed 2025
