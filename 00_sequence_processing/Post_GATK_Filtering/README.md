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
# Evaluate filtering using MSI's Open On Demand with code in this Jupyter Notebook
Evaluate_filtering-snps_introgression_project.ipynb
```

Additional filtering checks.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering
# Generates bcftools stats, calculate Ts/Tv, and calculates MAF
sbatch evaluate_vcf_filtering.sh
sbatch tstv_snps.sh
sbatch summarize_missingness.sh
```

Check het/hom_alt ratio for any outliers (generated from bcftools stats output file with some custom calculations in script).

```bash
# In dir: ~/Alignments/introgression_project/all_dom_and_wild/Filtered/vcf_summary
(head -n 1 dom_and_wild_snps_biallelic.callable.cap50x.het_hom_ratios.txt && tail -n +2 dom_and_wild_snps_biallelic.callable.cap50x.het_hom_ratios.txt | sort -k4,4nr) | head
PI_231151	36313	23903	0.658249
PI_564666	43936	19500	0.443827
PI_327716	34240	12854	0.375409
PI_327606	74009	19817	0.267765
PI_422230	74614	19488	0.261184
HOR14153	80747	16812	0.208206
BCC003	78043	12145	0.155619
WBDC_113	1561	157	0.100577
WBDC_095	1118	111	0.0992844
WBDC_115	1870	176	0.0941176
```

Check individuals with the most missing genotypes.

```bash
# In dir: ~/Alignments/introgression_project/all_dom_and_wild/Filtered/ann_visualization
sort -k5,5nr dom_and_wild_snps_biallelic.callable_missingness.imiss | head
```

Exclude the following individuals with very high proportion missing data (97%-99% missing).

```bash
WBDC_146
WBDC_095
WBDC_113
WBDC_115
WBDC_048
WBDC_122
```

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Post_GATK_Filtering
sbatch exclude_samples.sh
```
