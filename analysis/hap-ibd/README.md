# hap-ibd

Run hap-ibd (https://github.com/browning-lab/hap-ibd) to get a first pass guess at which wild and domesticated individuals share large segments of HBD/IBD segments. This can help us decide which wild samples to include in our query sample set when running FLARE for local ancestry inference.

### Prepare files

Combine all domesticated and wild (including likely introgressed wild) samples into a single VCF.

```bash
module load bcftools/1.10.2
module load htslib/1.9

# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o dom_and_wild_with_introgressed_merged.phased.imputed.vcf.gz
tabix -p vcf --csi dom_and_wild_with_introgressed_merged.phased.imputed.vcf.gz

# Exclude SNPs with missing genotypes
bcftools view -e "GT='mis'" dom_and_wild_with_introgressed_merged.phased.imputed.vcf.gz -O z -o dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz
tabix -p vcf --csi dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz
```

### Run hap-ibd

```bash
./run_hap-ibd-wbdc_nsgc_geno.sh
```

Pull out only pairs with WBDC samples since that's what we're primarily interested in.

```bash
# In dir: ~/Projects/Introgressed/hap_ibd
zcat dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz | grep "WBDC" | grep "PI"
```

Did not seem to make domesticated vs wild comparisons.
