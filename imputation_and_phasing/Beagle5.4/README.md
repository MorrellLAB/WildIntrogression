# Phasing and imputation with Beagle 5.4

Phasing and imputation with [Beagle 5.4](http://faculty.washington.edu/browning/beagle/beagle.html).

## Prepare input files

Files needed:

#### Separate domesticated samples and wild samples before phasing and imputation. Here are the filepaths:

```bash
# NSGC Core accessions
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/domesticated_snps.polymorphic.filt_miss_het.vcf.gz
# Wild accessions
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz
```

#### Plink MAP file

Prepare Plink MAP files for domesticated and wild vcfs. The best we can do is extrapolation from the BOPA and 9k consensus maps (https://github.com/MorrellLAB/morex_reference/tree/master/genetic_maps).

Convert VCF to Plink files and prepare Plink map files.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
# Dependencies
module load plink/1.90b6.10
module load python3/3.8.3_anaconda2020.07_mamba

# Convert VCF to Plink files and add cM positions to map files
# However, using plink to add cM positions doesn't account for SNPs without a cM position
#plink --vcf wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz --allow-extra-chr --keep-allele-order --update-cm ~/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt 3 2 --recode --make-bed --out wbdc_bopa_snps.polymorphic.filt_miss_het
plink --vcf wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz --allow-extra-chr --keep-allele-order --recode --make-bed --out wbdc_bopa_snps.polymorphic.filt_miss_het

# Tried interpolation script: https://github.com/MorrellLAB/Barley_Inversions/blob/master/01_analyses/Alchemy_genotype/script/Interpolate_MAP_Positions.py
# But didn't end up working, still a lot of missing cM positions and didn't get more SNPs to work compared
# to what we had

# Sort by physical locations in preparation for interpolating cM locations
#sort -k1,1 -k4,4n wbdc_bopa_snps.polymorphic.filt_miss_het.map > wbdc_bopa_snps.polymorphic.filt_miss_het.sorted_by_phys.map

# Interpolate missing cM locations

# Check for duplicate markers
plink --bfile wbdc_bopa_snps.polymorphic.filt_miss_het --allow-extra-chr --keep-allele-order -list-duplicate-vars ids-only suppress-first --out wbdc_duplicate_vars

# Add cM positions to Plink .map file
~/GitHub/WildIntrogression/imputation_and_phasing/Beagle5.4/add_cM_to_MAP_exclude_missing.py wbdc_bopa_snps.polymorphic.filt_miss_het.map ~/GitHub/morex_reference/genetic_maps/BOPA/2013_iSelect_Genetic_Map.txt 3 2 | sort -k1,1 -k3,3n -k4,4n > wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.map

# Identify SNPs where cM order doesn't agree with physical position ascending order
# (these cause errors for Beagle)
sort -k1,1 -k4,4n wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.map > wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.sort_by_physpos.map
# Find rows where sort order differs
diff wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.map wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.sort_by_physpos.map | grep '>\|<' | cut -f 2 | sort -uV > discordant_marker_order-wbdc.txt
# Remove likely problematic markers for Beagle
grep -vf discordant_marker_order-wbdc.txt wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.map | grep -vf wbdc_duplicate_vars.dupvar > wbdc_bopa_snps.polymorphic.filt_miss_het.updated_cm.excluded_problem_markers.map
```

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
# Convert VCF to Plink files and add cM positions to map files
plink --vcf domesticated_snps.polymorphic.filt_miss_het.vcf.gz --allow-extra-chr --keep-allele-order --recode --make-bed --out domesticated_snps.polymorphic.filt_miss_het

# Check for duplicate markers
plink --bfile domesticated_snps.polymorphic.filt_miss_het --allow-extra-chr --keep-allele-order -list-duplicate-vars ids-only suppress-first --out domesticated_duplicate_vars

# Add cM positions to Plink .map file
~/GitHub/WildIntrogression/imputation_and_phasing/Beagle5.4/add_cM_to_MAP_exclude_missing.py domesticated_snps.polymorphic.filt_miss_het.map ~/GitHub/morex_reference/genetic_maps/9k_iSelect_Genetic_Map.txt 4 1 | sort -k1,1 -k3,3n -k4,4n > domesticated_snps.polymorphic.filt_miss_het.updated_cm.map

# Identify SNPs where cM order doesn't agree with physical position ascending order
# (these cause errors for Beagle)
sort -k1,1 -k4,4n domesticated_snps.polymorphic.filt_miss_het.updated_cm.map > domesticated_snps.polymorphic.filt_miss_het.updated_cm.sort_by_physpos.map
# Find rows where sort order differs
diff domesticated_snps.polymorphic.filt_miss_het.updated_cm.map domesticated_snps.polymorphic.filt_miss_het.updated_cm.sort_by_physpos.map | grep '>\|<' | cut -f 2 | sort -uV > discordant_marker_order-domesticated.txt
# Remove likely problematic markers for Beagle
grep -vf discordant_marker_order-domesticated.txt domesticated_snps.polymorphic.filt_miss_het.updated_cm.map | grep -vf domesticated_duplicate_vars.dupvar > domesticated_snps.polymorphic.filt_miss_het.updated_cm.excluded_problem_markers.map
```

## Beagle phasing and imputation

```bash
# In dir: ~/GitHub/WildIntrogression/imputation_and_phasing/Beagle5.4
./run_beagle-wbdc_geno.sh
./run_beagle-nsgc_geno.sh
```

Check het after imputation

```bash
./check_het-wbdc_imputed.sh
./check_het-nsgc_imputed.sh
```
