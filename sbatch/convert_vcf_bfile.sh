#!/bin/bash
#SBATCH --job-name=convert_vcf_bfile
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/convert_vcf_bfile_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/convert_vcf_bfile_%A.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

i=${ARGS1}

plink2 --vcf /scratch/jll1/Bangladesh/liftOver/chr${i}.vcf --chr ${i} --make-bed --out /scratch/jll1/Bangladesh/bfile/raw_liftOver/chr${i}
