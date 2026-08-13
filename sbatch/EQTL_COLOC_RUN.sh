#!/bin/bash
#SBATCH --job-name=EQTL_COLOC_RUN
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/EQTL_COLOC_RUN_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/EQTL_COLOC_RUN_%A.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env

# importing imqtl extract list
extract_list=${ARGS1}

# running code to perform eqtl-imqtl coloc analysis
Rscript /gpfs/data/pierce-lab/james.li/imQTL/code/rscripts/EQTL_COLOC_RUN.R ${extract_list}
Rscript /gpfs/data/pierce-lab/james.li/imQTL/code/rscripts/EQTL_COLOC_RUN_custom_prior.R ${extract_list}