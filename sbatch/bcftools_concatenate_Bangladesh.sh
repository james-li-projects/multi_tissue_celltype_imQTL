#!/bin/bash
#SBATCH --job-name=bcftools_merge_ALL_DATASET
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/bcftools_concatenate_Bangladesh_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/bcftools_concatenate_Bangladesh_%A.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=250gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load bcftools/1.17

# setting number of processes 
ulimit -n 10000
# printing out this set number of processes
ulimit -n 

i=${ARGS1}

export TMPDIR=/scratch/jll1/Bangladesh/tmp/chr${i}
in_path=/scratch/jll1/Bangladesh

# indexing this merged/processed vcf
bcftools index -t ${out_path}/unfilt_chr${i}.vcf.gz


# concatenating the chromosomes
bcftools concat \
    ${in_path}/AABCG_dataset_01_WGS_SNP_chr${i}.multi.norm.vcf.gz \
    ${in_path}/AABCG_dataset_01_WGS_indel_chr${i}.multi.norm.vcf.gz | bcftools sort -Oz -o ${out_path}/AABCG_dataset_01_WGS_chr${i}.vcf.gz

# indexing the concatenated vcfs for each chromosome
bcftools index -t ${out_path}/AABCG_dataset_01_WGS_chr${i}.vcf.gz



