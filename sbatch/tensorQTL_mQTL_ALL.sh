#!/bin/bash
#SBATCH --job-name=tensorQTL_mQTL_BEST_ONLY
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/tensorQTL_mQTL_ALL_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/tensorQTL_mQTL_ALL_%A.err
#SBATCH --time=24:00:00
#SBATCH -p gpuq
#SBATCH --gres gpu:1
#SBATCH --mem-per-gpu=100gb

# loading modules 
module load gcc/12.1.0
module load python/3.10.5

# importing tissue and cell type arguments
dataset=${ARGS1}
tissue=${ARGS2}

# genetic data directory
genetic_data_dir=/gpfs/data/pierce-lab/james.li/imQTL/data/${dataset}/genetic_data

# covariate/methylation input directory
input_dir=/gpfs/data/pierce-lab/james.li/imQTL/input/${dataset}

# identifying input files and paths
plink_file=${genetic_data_dir}/processed_genetic_data_chrprefix
DNAm_bed=${input_dir}/${tissue}.bed
covariates_file=${input_dir}/processed_covariates_${tissue}.txt

# setting output file prefix
prefix="/gpfs/data/pierce-lab/james.li/imQTL/output/${dataset}/mQTL/full_assoc/tensorQTL_mQTL_${tissue}"

# mapping imQTLs
python3 -m tensorqtl ${plink_file} ${DNAm_bed} ${prefix} \
    --covariates ${covariates_file} \
    --window 500000 \
    --mode cis_nominal \
    --maf_threshold 0.1 \
    --seed 2024
