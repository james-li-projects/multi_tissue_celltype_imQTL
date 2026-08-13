#!/bin/bash
#SBATCH --job-name=GWAS_COLOC_EXTRACT_MQTL_DATA
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/GWAS_COLOC_EXTRACT_MQTL_DATA_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/GWAS_COLOC_EXTRACT_MQTL_DATA_%A.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env

# importing imqtl extract list
extract_list=${ARGS1}

# running code to plot qqplots
Rscript /gpfs/data/pierce-lab/james.li/imQTL/code/rscripts/GWAS_COLOC_EXTRACT_MQTL_DATA.R ${extract_list}
