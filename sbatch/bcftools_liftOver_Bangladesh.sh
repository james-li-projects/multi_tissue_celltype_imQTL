#!/bin/bash
#SBATCH --job-name=bcftools_liftOver_Bangladesh
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/bcftools_liftOver_Bangladesh_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/bcftools_liftOver_Bangladesh_%A.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load miniconda3
source activate pytorch_env
module load bcftools

# setting number of processes 
ulimit -n 10000
# printing out this set number of processes
ulimit -n 

i=${ARGS1}

export TMPDIR=/scratch/jll1/Bangladesh/tmp/chr${i}
in_path=/scratch/jll1/Bangladesh
out_path=/scratch/jll1/Bangladesh/bcftools_liftOver

# running bcftools
cd /scratch/jll1/imQTL/genetic_data/
hg19=/gpfs/data/pierce-lab/james.li/hg19/hg19.fa
hg38=/gpfs/data/pierce-lab/james.li/hg38/hg38.fa

input_chain=/gpfs/data/pierce-lab/james.li/liftOver/hg19ToHg38.over.chain.gz
input_vcf=${in_path}/unfilt_chr${i}.vcf.gz


bcftools +liftover ${input_vcf} -s ${hg19} -f ${hg38} -c ${input_chain} | bcftools sort -Ov -o ${out_path}/chr${i}.vcf
