#!/bin/bash
#SBATCH --job-name=bcftools_merge_ALL_DATASET
#SBATCH --output=/gpfs/data/pierce-lab/james.li/imQTL/logs/CrossMap_liftOver_Bangladesh_%A.out
#SBATCH --error=/gpfs/data/pierce-lab/james.li/imQTL/logs/CrossMap_liftOver_Bangladesh_%A.err
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
out_path=/scratch/jll1/Bangladesh/liftOver
fasta_file=/gpfs/data/pierce-lab/james.li/hg38/hg38.fa
input_chain=/gpfs/data/pierce-lab/james.li/liftOver/hg19ToHg38.over.chain.gz
unfiltered_input_vcf=${in_path}/unfilt_chr${i}.vcf.gz
renamedchr_input_vcf=${in_path}/renamedchr_chr${i}.vcf.gz

# renaming chromosomes
rename_chr_file=/gpfs/data/pierce-lab/james.li/bcftools_chr_name_conv/chr_name_conv.txt
bcftools annotate --rename-chrs ${rename_chr_file} ${unfiltered_input_vcf} -Oz -o ${renamedchr_input_vcf}

# running CrossMap
cd /scratch/jll1/imQTL/genetic_data/
CrossMap vcf ${input_chain} \
  ${renamedchr_input_vcf} \
  ${fasta_file} \
  ${out_path}/chr${i}.vcf
