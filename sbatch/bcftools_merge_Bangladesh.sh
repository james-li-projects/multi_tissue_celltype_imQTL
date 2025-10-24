#!/bin/bash
#SBATCH --job-name=bcftools_merge_ALL_DATASET
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/bcftools_merge_Bangladesh_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/bcftools_merge_Bangladesh_%A.err
#SBATCH --time=72:00:00
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
out_path=/scratch/jll1/Bangladesh

in_path=/gpfs/data/pierce-lab/james.li/imQTL/Bangladesh

bcftools index -t ${in_path}/custom/chr${i}.dose.vcf.gz
bcftools index -t ${in_path}/prev_w_exom/chr${i}.dose.vcf.gz
bcftools index -t ${in_path}/prev_wo_exom/chr${i}.dose.vcf.gz

# merging vcfs, creating specific variant IDs, removing duplicate variants 
bcftools merge -m none \
    ${in_path}/custom/chr${i}.dose.vcf.gz \
    ${in_path}/prev_w_exom/chr${i}.dose.vcf.gz \
    ${in_path}/prev_wo_exom/chr${i}.dose.vcf.gz -Oz -o ${out_path}/unfilt_chr${i}.vcf.gz

# indexing this merged/processed vcf
bcftools index -t ${out_path}/unfilt_chr${i}.vcf.gz
