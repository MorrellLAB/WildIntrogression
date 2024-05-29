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

#### Separate wild vs wild_introgressed (these are samples we suspect are introgressed).

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

---

### Run FLARE with additional potentially introgressed WBDC accessions

Update wild samples from `~/GitHub/WildIntrogression/analysis/hap-ibd/potentially_introgressed_wbdc_list.minmarkers40_maxgap1000.txt` likely introgressed vs not likely introgressed in primary spreadsheet `~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx` the tabs `wild_pop_info` and `all_combined_v2`, then pull sample names.

Modify the file: `wild_all_samples_groupings`

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/FLARE
# Generate sample name lists for introgressed vs not introgressed
grep "wild_introgressed" wild_all_samples_groupings.txt | cut -f 1 | sort -V > wild_likely_introgressed_samples.txt

grep -v "wild_introgressed" wild_all_samples_groupings.txt | cut -f 1 | sort -V > wild_samples_not_introgressed.txt
```

Use updated sample lists and separate wild from wild_introgressed in VCF.

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
48

bcftools query -l wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz | wc -l
270
```

Merge domesticated and wild (not introgressed) to use as reference panel VCF.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.vcf.gz

# Exclude SNPs with missing genotypes
bcftools view -e "GT='mis'" dom_and_wild_merged.phased.imputed.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.no_missing.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.no_missing.vcf.gz

# Check number of SNPs
zgrep -v "#" dom_and_wild_merged.phased.imputed.no_missing.vcf.gz | wc -l
1749
```

Re-run FLARE

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/FLARE
./run_flare-wbdc_geno.sh
```

Plot output with script:

```bash
plot_flare_output.R
```

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare
bedtools intersect -header -a ~/Dropbox/Projects/Wild_Introgression/Data/morex_v3_gff/Hv_Morex.pgsb.Jul2020.HC.gff3 -b plots/wbdc_likely_introgressed_segments-breeding.bed | uniq > Hv_Morex.pgsb.Jul2020.HC.wbdc_likely_introgressed_segments.gff3
```

---

### Run FLARE with three subsampled sets of 48 wild samples each as query

Purpose: comparison of introgressed segments for wild vs wild-introgressed samples

Randomly sample 48 wild samples (exclude wild-introgressed from this analysis)

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/FLARE
gshuf -n 48 wild_samples_not_introgressed.txt | sort -V > wild_not_introgressed-rand1_48_samp.txt
gshuf -n 48 wild_samples_not_introgressed.txt | sort -V > wild_not_introgressed-rand2_48_samp.txt
gshuf -n 48 wild_samples_not_introgressed.txt | sort -V > wild_not_introgressed-rand3_48_samp.txt
cp wild_not_introgressed-rand1_48_samp.txt ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
cp wild_not_introgressed-rand2_48_samp.txt ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
cp wild_not_introgressed-rand3_48_samp.txt ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
```

Merge domesticated and wild (not introgressed) to use as reference panel VCF, except in this analysis remove the randomly sampled wild to be used as the query samples.

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
# Set 1
bcftools view --samples-file wild_not_introgressed-rand1_48_samp.txt wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_not_introgressed-rand1_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.rand1_48_samp.vcf.gz
bcftools view --samples-file ^wild_not_introgressed-rand1_48_samp.txt wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_not_introgressed-rand1_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set1.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand1_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.rand1_48_samp.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand1_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set1.vcf.gz

# Set 2
bcftools view --samples-file wild_not_introgressed-rand2_48_samp.txt wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_not_introgressed-rand2_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.rand2_48_samp.vcf.gz
bcftools view --samples-file ^wild_not_introgressed-rand2_48_samp.txt wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_not_introgressed-rand2_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set2.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand2_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.rand2_48_samp.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand2_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set2.vcf.gz

# Set 3
bcftools view --samples-file wild_not_introgressed-rand3_48_samp.txt wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_not_introgressed-rand3_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.rand3_48_samp.vcf.gz
bcftools view --samples-file ^wild_not_introgressed-rand3_48_samp.txt wild.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o wild_not_introgressed-rand3_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set3.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand3_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.rand3_48_samp.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand3_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set3.vcf.gz
```

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
# Set 1
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wild_not_introgressed-rand1_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set1.vcf.gz -O z -o wild_not_introgressed-rand1_48_samp/dom_and_wild_merged.phased.imputed.set1.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand1_48_samp/dom_and_wild_merged.phased.imputed.set1.vcf.gz

# Exclude SNPs with missing genotypes
cd wild_not_introgressed-rand1_48_samp
bcftools view -e "GT='mis'" dom_and_wild_merged.phased.imputed.set1.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.no_missing.set1.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.no_missing.set1.vcf.gz

# Check number of SNPs
zgrep -v "#" dom_and_wild_merged.phased.imputed.no_missing.set1.vcf.gz | wc -l
    1749
```

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
# Set 2
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wild_not_introgressed-rand2_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set2.vcf.gz -O z -o wild_not_introgressed-rand2_48_samp/dom_and_wild_merged.phased.imputed.set2.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand2_48_samp/dom_and_wild_merged.phased.imputed.set2.vcf.gz

# Exclude SNPs with missing genotypes
cd wild_not_introgressed-rand2_48_samp
bcftools view -e "GT='mis'" dom_and_wild_merged.phased.imputed.set2.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.no_missing.set2.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.no_missing.set2.vcf.gz

# Check number of SNPs
zgrep -v "#" dom_and_wild_merged.phased.imputed.no_missing.set2.vcf.gz | wc -l
    1749
```

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
# Set 3
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wild_not_introgressed-rand3_48_samp/wild.polymorphic.filt_miss_het.phased.imputed.set3.vcf.gz -O z -o wild_not_introgressed-rand3_48_samp/dom_and_wild_merged.phased.imputed.set3.vcf.gz
tabix -p vcf --csi wild_not_introgressed-rand3_48_samp/dom_and_wild_merged.phased.imputed.set3.vcf.gz

# Exclude SNPs with missing genotypes
cd wild_not_introgressed-rand3_48_samp
bcftools view -e "GT='mis'" dom_and_wild_merged.phased.imputed.set3.vcf.gz -O z -o dom_and_wild_merged.phased.imputed.no_missing.set3.vcf.gz
tabix -p vcf --csi dom_and_wild_merged.phased.imputed.no_missing.set3.vcf.gz

# Check number of SNPs
zgrep -v "#" dom_and_wild_merged.phased.imputed.no_missing.set3.vcf.gz | wc -l
    1749
```

Prepare reference panel

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/flare_input_files
# For this analysis, exlude WBDC samples identified as wild-introgressed by FLARE
grep "wild" dom_and_wild_not_introgressed.ref.panel | grep -vf ~/GitHub/WildIntrogression/analysis/FLARE/wild_likely_introgressed_samples.txt > wild_not_introgressed.ref.panel
# Set 1
grep -v "wild" dom_and_wild_not_introgressed.ref.panel > wild_not_introgressed-rand1_48_samp/dom_and_wild_not_introgressed.rand1_48samp.ref.panel
grep -vf wild_not_introgressed-rand1_48_samp.txt wild_not_introgressed.ref.panel >> wild_not_introgressed-rand1_48_samp/dom_and_wild_not_introgressed.rand1_48samp.ref.panel

# Set 2
grep -v "wild" dom_and_wild_not_introgressed.ref.panel > wild_not_introgressed-rand2_48_samp/dom_and_wild_not_introgressed.rand2_48samp.ref.panel
grep -vf wild_not_introgressed-rand2_48_samp.txt wild_not_introgressed.ref.panel >> wild_not_introgressed-rand2_48_samp/dom_and_wild_not_introgressed.rand2_48samp.ref.panel

# Set 3
grep -v "wild" dom_and_wild_not_introgressed.ref.panel > wild_not_introgressed-rand3_48_samp/dom_and_wild_not_introgressed.rand3_48samp.ref.panel
grep -vf wild_not_introgressed-rand3_48_samp.txt wild_not_introgressed.ref.panel >> wild_not_introgressed-rand3_48_samp/dom_and_wild_not_introgressed.rand3_48samp.ref.panel
```

Run FLARE

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/FLARE
# Set 1
./run_flare-wbdc_geno_rand1_48samp.sh
# Copy results to Dropbox - had errors calling on files stored in Dropbox
cp -r ~/Downloads/temp_flare/wild_not_introgressed-rand1_48_samp_out ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare

# Set 2
./run_flare-wbdc_geno_rand2_48samp.sh
cp -r ~/Downloads/temp_flare/wild_not_introgressed-rand2_48_samp_out ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare

# Set 3
./run_flare-wbdc_geno_rand3_48samp.sh
cp -r ~/Downloads/temp_flare/wild_not_introgressed-rand3_48_samp_out ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare
```
