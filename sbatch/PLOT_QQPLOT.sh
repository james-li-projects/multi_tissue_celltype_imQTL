#!/bin/bash
#SBATCH --job-name=PLOT_QQPLOT
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/PLOT_QQPLOT_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/PLOT_QQPLOT_%A.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate pytorch_env

# importing tissue and cell type arguments
dataset=${ARGS1}
combination=${ARGS2}

# running code to plot qqplots
python /gpfs/data/pierce-lab/james.li/imQTL/code/python/PLOT_QQPLOT.py ${dataset} ${combination}
