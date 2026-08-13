#!/bin/bash
#SBATCH --job-name=GWAS_COLOC_BULKMQTL_EXTRACT_WINDOWS
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/GWAS_COLOC_BULKMQTL_EXTRACT_WINDOWS_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/GWAS_COLOC_BULKMQTL_EXTRACT_WINDOWS_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env

# importing argument for GWAS file
mqtl_chunk=${ARGS1}

# running code to perform imqtl-gwas coloc analysis
Rscript /gpfs/data/pierce-lab/james.li/imQTL/code/rscripts/GWAS_COLOC_BULKMQTL_EXTRACT_WINDOWS.R ${mqtl_chunk}
