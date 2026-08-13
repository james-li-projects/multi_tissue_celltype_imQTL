#!/bin/bash
#SBATCH --job-name=OBTAIN_MASHR_INPUTS
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/OBTAIN_MASHR_INPUTS_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/OBTAIN_MASHR_INPUTS_%A.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate pytorch_env

# importing tissue and cell type arguments
dataset=${ARGS1}
combination=${ARGS2}

# running code to plot qqplots
python /gpfs/data/pierce-lab/james.li/imQTL/code/python/OBTAIN_MASHR_INPUTS.py ${dataset} ${combination}
