# PCA Analsis

Do a PCA analysis for the WBDC and domesticated lines genotype data.

### Plink PCA

```bash
# In dir: ~/Projects/Introgressed/plink_pca
# Dependencies
module load plink/1.90b6.10

vcf_file="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.final.vcf.gz"
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/plink_pca"
out_prefix="dom_and_wbdc_exome_capture"

# Run plink PCA
plink --vcf ${vcf_file} --pca --keep-allele-order --allow-extra-chr --double-id --out ${out_dir}/${out_prefix}.pca
```

Run Plink PCA on just the wild barley genotyped samples.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
vcf_file="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz"
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/plink_pca"
out_prefix="wbdc_geno_samples"

# Run plink PCA
plink --vcf ${vcf_file} --pca --keep-allele-order --allow-extra-chr --out ${out_dir}/${out_prefix}.pca
```

Use `plot_plink_PCA.R` script to generate plots.


Run Plink PCA on the wild and domesticated genotyped samples.

```bash
vcf_file="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz"
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/plink_pca"
out_prefix="dom_and_wbdc_geno_samples"

# Run plink PCA
plink --vcf ${vcf_file} --pca --keep-allele-order --allow-extra-chr --out ${out_dir}/${out_prefix}.pca
```

Run Plink PCA on the WBDC GBS dataset.

```bash
vcf_file="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/Ahmad_GBS_morex_v3/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz"
out_dir="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/plink_pca"
out_prefix="wbdc_gbs"

# Run plink PCA
plink --vcf ${vcf_file} --pca --keep-allele-order --allow-extra-chr --out ${out_dir}/${out_prefix}.pca
```

Plot PCA using final wild introgressed assignments listed in Table S1 based on all analyses (including FLARE, ABBABABA, etc.).

```bash
plot_plink_PCA_flare_assign.R
```
