#!/bin/bash
#SBATCH --job-name=PLINK2_PCA
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/PLINK2_PCA.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/PLINK2_PCA.err

module load gcc/12.1.0
module load plink/2.0

DATASET=${ARGS1}

# pruning dataset
plink2 --bfile /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/processed_genetic_data \
  --indep-pairwise 500kb 0.2 \
  --threads 1 \
  --memory 100000 \
  --out /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/pruned

# obtaining a pfile that has been pruned
plink2 --bfile /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/processed_genetic_data \
  --extract /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/pruned.prune.in \
  --threads 1 \
  --memory 100000 \
  --make-bed \
  --out /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/PRUNED

# generating PCs
plink2 --bfile /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/PRUNED \
    --pca 30 \
    --threads 16 \
    --memory 100000 \
    --out /gpfs/data/pierce-lab/james.li/imQTL/data/${DATASET}/genetic_data/PLINK2_PCA_30
    
