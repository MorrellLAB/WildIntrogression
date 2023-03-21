# PCA Analsis

Do a PCA analysis for the WBDC and domesticated lines genotype data.

### Plink PCA

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_pca
vcf_file="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf"
out_dir="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_pca"
out_prefix="dom_and_wbdc_samples"

# Dependencies
module load plink/1.90b6.10

# Run plink PCA
plink --vcf ${vcf_file} --pca --allow-extra-chr --out ${out_dir}/${out_prefix}.pca
```

Run Plink PCA on just the wild barley samples.

```bash
# In dir: /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_pca
vcf_file="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/wbdc_318_BOPA_morex_v3.vcf"
out_dir="/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/plink_pca"
out_prefix="wbdc_samples"

# Dependencies
module load plink/1.90b6.10

# Run plink PCA
plink --vcf ${vcf_file} --pca --allow-extra-chr --out ${out_dir}/${out_prefix}.pca
```

Use `plot_plink_PCA.R` script to generate plots.
