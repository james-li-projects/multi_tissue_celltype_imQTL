
#####################################
# initializing chromatin segmentation states
chromatin_states=(15_Quies 1_TssA 14_ReprPCWk 9_Het 7_Enh 5_TxWk 2_TssAFlnk 10_TssBiv 4_Tx 13_ReprPC 12_EnhBiv 11_BivFlnk 8_ZNF/Rpts 6_EnhG 3_TxFlnk)

#####################################
# obtaining colon chromatin states in hg19
input_file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/download_data/chromatin_segmentation/E075_15_coreMarks_dense.bed.gz"
output_dir="/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg19"
for chromatin_state in "${chromatin_states[@]}"; do
  safe_name=$(echo "$chromatin_state" | sed 's|/|_|g')
  zgrep "$chromatin_state" "$input_file" > "${output_dir}/colon_${safe_name}.bed"
done

#####################################
# obtaining lung chromatin states in hg19
input_file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/download_data/chromatin_segmentation/E096_15_coreMarks_dense.bed.gz"
output_dir="/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg19"
for chromatin_state in "${chromatin_states[@]}"; do
  safe_name=$(echo "$chromatin_state" | sed 's|/|_|g')
  zgrep "$chromatin_state" "$input_file" > "${output_dir}/lung_${safe_name}.bed"
done

#####################################
# obtaining whole blood chromatin states in hg19
#####################################
# E039 Primary T helper naive cells #
# E038 Primary T helper naive cells #
# E047   Primary T CD8+ naive cells #
# E029            Primary monocytes #
# E032              Primary B cells #
# E046 Primary natural killer cells #
# E030          Primary neutrophils #
#####################################
rm /gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg19/wb_*.bed
for E_ID in E039 E038 E047 E029 E032 E046 E030
do
  echo Processing $E_ID
  input_file="/gpfs/data/phs/groups/Projects/GEMS/james.li/EWAS/manuscript/download_data/chromatin_segmentation/${E_ID}_15_coreMarks_dense.bed.gz"
output_dir="/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg19"
for chromatin_state in "${chromatin_states[@]}"; do
  safe_name=$(echo "$chromatin_state" | sed 's|/|_|g')
  zgrep "$chromatin_state" "$input_file" >> "${output_dir}/wb_${safe_name}.bed"
done
done

# converting bed files to hg38 using liftOver
liftOver_software="/gpfs/data/pierce-lab/james.li/imQTL/software/liftOver/liftOver"
UCSC_liftOver_chain="/gpfs/data/pierce-lab/james.li/imQTL/software/liftOver/hg19ToHg38.over.chain.gz"
output_dir_mapped=/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg38/mapped
output_dir_unmapped=/gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg38/unmapped
cd /gpfs/data/pierce-lab/james.li/imQTL/data/ADDITIONAL/CHROMATIN_SEGMENTATION/hg19
for input_file_prefix in *
do
  echo ${input_file_prefix}
  ${liftOver_software} ${input_file_prefix} ${UCSC_liftOver_chain} ${output_dir_mapped}/${input_file_prefix} ${output_dir_unmapped}/${input_file_prefix} 
done
