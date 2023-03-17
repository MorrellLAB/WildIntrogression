# Post GATK Filtering

We followed the GATK best practices pipeline for SNP calling implemented in [`sequence_handling`](https://github.com/MorrellLAB/sequence_handling) all the way through the handler `Variant_Recalibrator`. The post GATK filtering steps are documented below.

### Methods: SNPs

Keep only polymorphic sites.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering
sbatch Filter_VCF_Poly.sh
```

Next, filter on proportion heterozygous genotypes, min DP, proportion missing, quality, GQ, and filter by allelic balance. Also filter to biallelic sites.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering
sbatch Filter_VCF_Het_AB_DP_Miss_Qual_GQ.sh
```

### Check filtering

Add custom annotations for "known" vs "novel" variants and for each level of filtering ("filtered" vs "retained_filt1" vs "retained_filt2"). Then convert VCF to HDF5 format (the HDF5 format works better with scikit-allel).

- For SNPs a VCF of known variants (i.e., BOPA, 9K, 50K SNPs) was created. SNPs in the known variants VCF has the annotation "known", the remaining SNPs have the annotation "novel".
- There are too many SNPs to visualize all of them at once, so we will pull out chr1H only for visualization and do some downsampling. The purpose of visualization is to help tune filtering cutoffs and get a sense of how good of a job we are doing on filtering, so we don't need to include the full set of SNPs. One way to pull out just chr1H is using something like `bcftools view -t "chr1H_part1,chr1H_part2" raw_variants_snps.vcf > chr1H_raw_variants_snps.vcf`.

Prepare files for filtering evaluation.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering
# Add custom annotations
sbatch prep_ann_snps.sh
# Prepare chr1H VCF for visualization in Jupyter Notebooks
sbatch scikit_allel_plot_prep.sh
```
