#!/bin/bash
#SBATCH --job-name=GWAS_MQTL_COLOC_RUN
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/GWAS_MQTL_COLOC_RUN_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/GWAS_MQTL_COLOC_RUN_%A.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env

# importing argument for GWAS file
gwas_file=${ARGS1}

# running code to perform imqtl-gwas coloc analysis
Rscript /gpfs/data/pierce-lab/james.li/imQTL/code/rscripts/GWAS_MQTL_COLOC_RUN.R ${gwas_file}
