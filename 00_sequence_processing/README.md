# Data Processing

This directory contains config files and subdirectories for processing the exome capture resequencing data of wild and cultivated samples.

1. Exome capture resequencing data: config files were used to run sequence_handling. SNP filtering scripts are in subdirectory `Post_GATK_Filtering`
2. BOPA and 9K genotyping data: scripts and descriptions are in subdirectory `genotype_data_processing`
3. GBS data: see subdirectory `Ahmad_gbs_processing`

---

### Calculate coverage

Calculate the depth of each BAM file with mosdepth.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing
sbatch --array=0-634 run_mosdepth.sh
```

Summarize average genome-wide coverage per sample.

```bash
# In dir: ~/Projects/Introgressed/coverage_exome
module load datamash_ML/1.3
grep -v "chrUn" PI_498431.mosdepth.summary.txt | grep -v "total" | grep -v "_region"
```
