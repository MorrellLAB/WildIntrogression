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
