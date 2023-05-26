# Local ancestry inference with FLARE

Inferring local ancestry using [FLARE](https://github.com/browning-lab/flare).

## Prepare input files

### Prepare VCF

Prepare reference and admixed study sample VCFs.

For FLARE, All genotypes must be phased and have no missing alleles. We imputed and phased data using Beagle, see `imputation_and_phasing/Beagle5.4` directory for scripts/documentation.

### Prepare genetic map file

Prepare genetic map file, needs to be in PLINK `.map` format.

We'll combine the Plink .map files we used during genotype imputation.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
cat wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.map domesticated_snps.polymorphic.filt_miss_het.updated_cm.map | sort -k1,1 -k3,3n -k4,4n | uniq > wbdc_and_dom_snps.polymorphic.filt_miss_het.updated_cm.map

# Combine already identified duplicate vars
cat wbdc_duplicate_vars.dupvar domesticated_duplicate_vars.dupvar | sort -V > wbdc_and_dom_duplicate_vars.dupvar

# Identify SNPs where cM order doesn't agree with physical position ascending order
# (these cause errors for FLARE)
sort -k1,1 -k4,4n wbdc_and_dom_snps.polymorphic.filt_miss_het.updated_cm.map > wbdc_and_dom_snps.polymorphic.filt_miss_het.updated_cm.sort_by_physpos.map
# Find rows where sort order differs
diff wbdc_and_dom_snps.polymorphic.filt_miss_het.updated_cm.map wbdc_and_dom_snps.polymorphic.filt_miss_het.updated_cm.sort_by_physpos.map | grep '>\|<' | cut -f 2 | sort -uV > discordant_marker_order-wbdc_and_dom.txt
# Remove likely problematic markers for FLARE
grep -vf discordant_marker_order-wbdc_and_dom.txt wbdc_and_dom_snps.polymorphic.filt_miss_het.updated_cm.map | grep -vf wbdc_and_dom_duplicate_vars.dupvar > wbdc_and_dom_snps.polymorphic.filt_miss_het.excluded_problem_markers.map
```

### Prepare `ref-panel`

Sample groupings are pulled from the spreadsheet: `~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx`

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed
module load bcftools/1.10.2
module load htslib/1.9

bcftools query -l domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz | sort -V > sample_names_domesticated_snps.polymorphic.filt_miss_het.phased.imputed.txt
```

Check unique groupings in ref-panel.

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/FLARE
cut -f 2 domesticated_ref_panel_pops.txt | sort -uV
breeding
cultivar
genetic
landrace
uncertain
```

Separate wild vs wild_introgressed (these are samples we suspect are introgressed).

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed
# Select likely introgressed samples only
bcftools view --samples-file ~/GitHub/WildIntrogression/analysis/FLARE/wild_likely_introgressed_samples.txt wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_likely_introgressed.polymorphic.filt_miss_het.phased.imputed.vcf.gz
tabix -p vcf --csi wild_likely_introgressed.polymorphic.filt_miss_het.phased.imputed.vcf.gz

# Select wild (not introgressed)
bcftools view --samples-file ~/GitHub/WildIntrogression/analysis/FLARE/wild_samples_not_introgressed.txt wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz
tabix -p vcf --csi wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz

# Check sample counts
bcftools query -l wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz | wc -l
318

bcftools query -l wild_likely_introgressed.polymorphic.filt_miss_het.phased.imputed.vcf.gz | wc -l
13

bcftools query -l wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz | wc -l
305
```

Merge domesticated and wild (not introgressed) to use as reference panel VCF.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.vcf.gz

# Exclude SNPs with missing genotypes
bcftools view -e "GT='mis'" dom_and_wild_merged.phased.imputed.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.no_missing.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.no_missing.vcf.gz
```

`dom_and_wild_merged.phased.imputed.no_missing.vcf.gz`: 1,749 SNPs

### Run FLARE

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/FLARE
./run_flare-wbdc_geno.sh
```
