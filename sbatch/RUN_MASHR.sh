#!/bin/bash
#SBATCH --job-name=RUN_MASHR
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/RUN_MASHR_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/RUN_MASHR_%A.err
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env

# running code to perform mashr analysis
Rscript /gpfs/data/pierce-lab/james.li/imQTL/code/rscripts/ANALYSIS/RUN_MASHR.R 
