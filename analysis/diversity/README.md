# Fst

Calculate Fst using:

Exome capture dataset, prepare windows

```bash
module load bedtools/2.29.2

cd ~/Projects/Introgressed/diversity
bedtools makewindows -b /panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/captured_50x_morex_v3_partsRef.bed -w 600 -s 200 > captured_50x_morex_v3_partsRef.win600_step200.bed
```

Convert VCF to .geno.gz format

```bash

```

```bash
module load python3/3.8.3_anaconda2020.07_mamba
module load htslib/1.9
module load bcftools/1.10.2
GGDIR="/panfs/jay/groups/9/morrellp/liux1299/Software/genomics_general"

VCF="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz"
OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/diversity"
POPS_FILE=~/GitHub/WildIntrogression/analysis/ABBA-BABA/dom-wild.pops.txt

WIN_SIZE="10000"
STEP_SIZE="1000"
MIN_GS="100"

out_prefix=$(basename $VCF .vcf.gz)

# # Convert VCF to .geno.gz format
# python ${GGDIR}/VCF_processing/parseVCF.py -i ${VCF} | bgzip > ${OUT_DIR}/${out_prefix}.geno.gz

# Generate pops file
bcftools query -l ${VCF} > ${OUT_DIR}/vcf_sample_names_dom_and_wild.txt
~/GitHub/WildIntrogression/analysis/ABBA-BABA/generate_pops_file.py ~/GitHub/WildIntrogression/analysis/ABBA-BABA/Sample_Info_Complete.csv ${OUT_DIR}/vcf_sample_names_dom_and_wild.txt > ${OUT_DIR}/dom_wild.pops.txt

# Calculate sliding window diversity, we're interested in fst
#python ${GGDIR}/popgenWindows.py -w ${WIN_SIZE} -s ${STEP_SIZE} -m ${MIN_GS} -g ${OUT_DIR}/${out_prefix}.geno.gz -o ${OUT_DIR}/${out_prefix}.csv.gz -f phased --addWindowID --threads 5 -p "domesticated" -p "wild" -p "wild_introgressed" --popsFile ${POPS_FILE}
cd ~/GitHub/WildIntrogression/analysis/diversity
sbatch windowed_diversity_estimates.sh
```
