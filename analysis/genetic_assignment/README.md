# Genetic assignment: STRUCTURE analysis

We'll need to do LD pruning before starting.

### Data preparation

First, generate list of sample names and the population they belong to. Make sure these names stay in the same order as in the VCF file.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3
#zgrep "#CHROM" merged_domesticated_and_wbdc_318_morex_v3.vcf.gz | cut -f 10- | tr '\t' '\n' > sample_names_merged_domesticated_and_wbdc_318_morex_v3.txt
zgrep "#CHROM" wbdc_318_BOPA_morex_v3.vcf.gz | cut -f 10- | tr '\t' '\n' > sample_names_wbdc_318_BOPA_morex_v3.txt

# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
# Download sample list to local computer for processing
#scp liux1299@mesabi.msi.umn.edu:~/Projects/Introgressed/vcf/morex_v3/sample_names_merged_domesticated_and_wbdc_318_morex_v3.txt .
scp liux1299@mesabi.msi.umn.edu:~/Projects/Introgressed/vcf/morex_v3/sample_names_wbdc_318_BOPA_morex_v3.txt .
```

Add population information.

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
#./add_pop_info.py sample_names_merged_domesticated_and_wbdc_318_morex_v3.txt accession_types_codes.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "all_combined" > pop_info_v1.txt
./add_pop_info.py sample_names_wbdc_318_BOPA_morex_v3.txt accession_types_codes_wild.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "wild_pop_info" > pop_info_wild_only.txt

# Get a list of sample names where we don't have population info yet
#grep "NA" pop_info_v1.txt | cut -d':' -f 2 > missing_pop_info_sample_names.txt

# Re-run once we have filled in missing info on accession type
#./add_pop_info.py sample_names_merged_domesticated_and_wbdc_318_morex_v3.txt accession_types_codes.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "all_combined_v2" > pop_info_v1.txt
```

Rename VCF.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3
module load bcftools/1.10.2
#bcftools reheader -s ~/GitHub/WildIntrogression/analysis/genetic_assignment/pop_info_v1.txt -o ~/Projects/Introgressed/genetic_assignment/data/merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.vcf merged_domesticated_and_wbdc_318_morex_v3.vcf
bcftools reheader -s ~/GitHub/WildIntrogression/analysis/genetic_assignment/pop_info_wild_only.txt -o ~/Projects/Introgressed/genetic_assignment/data/wbdc_318_BOPA_morex_v3_wPopInfo.vcf wbdc_318_BOPA_morex_v3.vcf
```

Convert VCF to Plink 1.9 MAP/PED files.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
module load plink/1.90b6.10
# Convert VCF to Plink MAP/PED
#plink --vcf merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.vcf --allow-extra-chr --id-delim ":" --recode --out /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/data/merged_domesticated_and_wbdc_318_morex_v3_wPopInfo
plink --vcf wbdc_318_BOPA_morex_v3_wPopInfo.vcf --allow-extra-chr --id-delim ":" --recode --out /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/genetic_assignment/data/wbdc_318_BOPA_morex_v3_wPopInfo
```

In the MAP file, make sure chromosomes are designated with their chromosome number only and don't have additional characters. Example, `chr1H` gets changed to `1`.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
# We'll use sed to do an in-place edit
for i in $(seq 1 7)
do
    #sed -i "s/chr${i}H/${i}/g" merged_domesticated_and_wbdc_318_morex_v3_wPopInfo.map
    sed -i "s/chr${i}H/${i}/g" wbdc_318_BOPA_morex_v3_wPopInfo.map
done
```

Convert PED/MAP files to Plink BED/BIM files.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
#plink --file merged_domesticated_and_wbdc_318_morex_v3_wPopInfo --make-bed --out merged_domesticated_and_wbdc_318_morex_v3_wPopInfo
plink --file wbdc_318_BOPA_morex_v3_wPopInfo --make-bed --out wbdc_318_BOPA_morex_v3_wPopInfo
```

LD pruning with Plink 1.90b6.10.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
#plink --bfile merged_domesticated_and_wbdc_318_morex_v3_wPopInfo --indep-pairwise 50 10 0.1
plink --bfile wbdc_318_BOPA_morex_v3_wPopInfo --indep-pairwise 50 10 0.1
# Output
PLINK v1.90b6.10 64-bit (17 Jun 2019)          www.cog-genomics.org/plink/1.9/
(C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to plink.log.
Options in effect:
  --bfile wbdc_318_BOPA_morex_v3_wPopInfo
  --indep-pairwise 50 10 0.1

64123 MB RAM detected; reserving 32061 MB for main workspace.
2963 variants loaded from .bim file.
318 people (0 males, 0 females, 318 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 318 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.944915.
2963 variants and 318 people pass filters and QC.
Note: No phenotypes present.
Pruned 140 variants from chromosome 1, leaving 210.
Pruned 222 variants from chromosome 2, leaving 266.
Pruned 189 variants from chromosome 3, leaving 277.
Pruned 162 variants from chromosome 4, leaving 198.
Pruned 224 variants from chromosome 5, leaving 326.
Pruned 149 variants from chromosome 6, leaving 196.
Pruned 164 variants from chromosome 7, leaving 240.
Pruning complete.  1250 of 2963 variants removed.
Marker lists written to plink.prune.in and plink.prune.out .
```

There are 2963 - 1250 = 1,713 SNPs remaining for Admixture/STRUCTURE analysis.

Generate BIM/BED files after LD pruning.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
#plink --bfile merged_domesticated_and_wbdc_318_morex_v3_wPopInfo --extract plink.prune.in --make-bed --out pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo
plink --bfile wbdc_318_BOPA_morex_v3_wPopInfo --extract plink.prune.in --make-bed --out pruned_wbdc_318_BOPA_morex_v3_wPopInfo
```

### Running Structure

Convert to Structure format.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
plink --bfile pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo --recode structure --out pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo
```

### Running Admixture

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
sbatch run_admixture.sh
```

Find the K based on the lowest standard error.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/Admixture
grep -h CV log*.out | sort -V
# Output
CV error (K=2): 0.93943
CV error (K=3): 0.89593
CV error (K=4): 0.86695
CV error (K=5): 0.83790
CV error (K=6): 0.82559
CV error (K=7): 0.80119
CV error (K=8): 0.78878
CV error (K=9): 0.77984
CV error (K=10): 0.77100
CV error (K=11): 0.76811
CV error (K=12): 0.75587
CV error (K=13): 0.74774
CV error (K=14): 0.74152
CV error (K=15): 0.73701
CV error (K=16): 0.72565
CV error (K=17): 0.72366
CV error (K=18): 0.71855
CV error (K=19): 0.71186
CV error (K=20): 0.70625

# Save to file for plotting
grep -h CV log*.out | sort -V | awk '{ print $3,$4 }' | sed -e 's,(,,' -e 's,),,' -e 's,:,,' -e 's,=, ,' > cross-validation_error_estimates.txt
```
