# Genetic assignment: STRUCTURE analysis

We'll need to do LD pruning before starting.

### Data preparation

First, generate list of sample names and the population they belong to. Make sure these names stay in the same order as in the VCF file.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3
zgrep "#CHROM" wbdc_318_BOPA_morex_v3.vcf.gz | cut -f 10- | tr '\t' '\n' > sample_names_wbdc_318_BOPA_morex_v3.txt

# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
# Download sample list to local computer for processing
scp liux1299@mesabi.msi.umn.edu:~/Projects/Introgressed/vcf/morex_v3/sample_names_wbdc_318_BOPA_morex_v3.txt .
```

Add population information.

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
./add_pop_info.py sample_names_wbdc_318_BOPA_morex_v3.txt accession_types_codes_wild.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "wild_pop_info" > pop_info_wild_only.txt

# Get a list of sample names where we don't have population info yet
#grep "NA" pop_info_v1.txt | cut -d':' -f 2 > missing_pop_info_sample_names.txt

# Re-run once we have filled in missing info on accession type
#./add_pop_info.py sample_names_merged_domesticated_and_wbdc_318_morex_v3.txt accession_types_codes.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "all_combined_v2" > pop_info_v1.txt
```

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
./prep_wbdc_geno.sh
```

---

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
sbatch run_admixture-wild_geno.sh
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
