# Genetic assignment: STRUCTURE analysis

We'll need to do LD pruning before starting.

### Data preparation

#### Prepare WBDC genotype data only

First, generate list of sample names and the population they belong to. Make sure these names stay in the same order as in the VCF file.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
zgrep "#CHROM" wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz | cut -f 10- | tr '\t' '\n' > sample_names_wbdc_318_BOPA_morex_v3.txt

# Local computer
cd ~/GitHub/WildIntrogression/analysis/genetic_assignment
# Download sample list to local computer for processing
scp liux1299@mesabi.msi.umn.edu:~/Projects/Introgressed/vcf/morex_v3/Filtered/sample_names_wbdc_318_BOPA_morex_v3.txt .

# Add population info
./add_pop_info.py sample_names_wbdc_318_BOPA_morex_v3.txt accession_types_codes_wild.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "wild_pop_info" > pop_info_wild_only.txt

# Get a list of sample names where we don't have population info yet
#grep "NA" pop_info_v1.txt | cut -d':' -f 2 > missing_pop_info_sample_names.txt
```

Prepare VCF (includes LD pruning).

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
./prep_wbdc_geno.sh
```

#### Prepare NSGC and WBDC genotype data combined

First, generate list of sample names and the population they belong to. Make sure these names stay in the same order as in the VCF file.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered
zgrep "#CHROM" nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz | cut -f 10- | tr '\t' '\n' > sample_names_nsgc_and_wbdc_morex_v3.txt

# Local computer
cd ~/GitHub/WildIntrogression/analysis/genetic_assignment
# Download sample list to local computer for processing
scp liux1299@mesabi.msi.umn.edu:~/Projects/Introgressed/vcf/morex_v3/Filtered/sample_names_nsgc_and_wbdc_morex_v3.txt .

./add_pop_info.py sample_names_nsgc_and_wbdc_morex_v3.txt accession_types_codes.txt ~/Dropbox/Projects/Wild_Introgression/Data/Compiled_Introgression_Sample_List.xlsx "all_combined_v2" > pop_info_nsgc_and_wild_geno.txt
```

Prepare VCF (includes LD pruning).

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
./prep_nsgc_and_wbdc_geno.sh
```

---

### Running Admixture

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/genetic_assignment
# Utilize Slurm job arrays across replicate runs and GNU parallel across K values
sbatch --array=0-9 run_admixture-wbdc_geno.sh
sbatch --array=0-9 run_admixture-nsgc_and_wild_geno.sh
```

### Prepare population files for plotting

Prepare ind2pop.txt files for pong.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data_morex_v3
module load python3/3.8.3_anaconda2020.07_mamba

# WBDC BOPA snps only
~/GitHub/WildIntrogression/analysis/genetic_assignment/fam_to_pong_pop.py wbdc_bopa_snps_morex_v3_wPopInfo.pruned.fam ~/GitHub/WildIntrogression/analysis/genetic_assignment/accession_types_codes_wild.txt > ~/GitHub/WildIntrogression/analysis/genetic_assignment/ind2pop-wbdc_bopa_snps.txt

# NSGC and WBDC genotype data
~/GitHub/WildIntrogression/analysis/genetic_assignment/fam_to_pong_pop.py nsgc_and_wbdc_geno_snps_morex_v3_wPopInfo.pruned.fam ~/GitHub/WildIntrogression/analysis/genetic_assignment/accession_types_codes.txt > ~/GitHub/WildIntrogression/analysis/genetic_assignment/ind2pop-nsgc_and_wbdc_geno.txt
```

Prepare filemap for pong.

```bash
cd ~/Projects/Introgressed/genetic_assignment/Admixture_nsgc_and_wild
# Column3 paths must be relative to the directory of the .Qfilemap
find -name "*.Q" | sed 's,./,,' | tr '/' '\t' | tr '.' '\t' | awk '{ print "k"$4 $1 "\t" $4 "\t" $1"/"$2"."$3"."$4"."$5 }' | sort -k2,2V -k1,1V > nsgc_and_wbdc_geno_snps_morex_v3_wPopInfo.pruned.Qfilemap
head -n 90 nsgc_and_wbdc_geno_snps_morex_v3_wPopInfo.pruned.Qfilemap > nsgc_and_wbdc_geno_snps_morex_v3_wPopInfo.pruned.k2-10.Qfilemap

cd ~/Projects/Introgressed/genetic_assignment/Admixture_wild_only_pruned
# Column3 paths must be relative to the directory of the .Qfilemap
find -name "*.Q" | sed 's,./,,' | tr '/' '\t' | tr '.' '\t' | awk '{ print "k"$4 $1 "\t" $4 "\t" $1"/"$2"."$3"."$4"."$5 }' | sort -k2,2V -k1,1V > wbdc_bopa_snps_morex_v3_wPopInfo.pruned.Qfilemap
head -n 90 wbdc_bopa_snps_morex_v3_wPopInfo.pruned.Qfilemap > wbdc_bopa_snps_morex_v3_wPopInfo.pruned.k2-10.Qfilemap
```

### Run pong

On local computer, run pong (so we can open visualization results in browser):

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/genetic_assignment/Admixture_nsgc_and_wild
#pong -m nsgc_and_wbdc_geno_snps_morex_v3_wPopInfo.pruned.Qfilemap -i ~/GitHub/WildIntrogression/analysis/genetic_assignment/ind2pop-nsgc_and_wbdc_geno.txt
pong -m nsgc_and_wbdc_geno_snps_morex_v3_wPopInfo.pruned.k2-10.Qfilemap -i ~/GitHub/WildIntrogression/analysis/genetic_assignment/ind2pop-nsgc_and_wbdc_geno.txt

pong -m wbdc_bopa_snps_morex_v3_wPopInfo.pruned.k2-10.Qfilemap -i ~/GitHub/WildIntrogression/analysis/genetic_assignment/ind2pop-wbdc_bopa_snps.txt
```

---
### Running Structure

Convert to Structure format.

```bash
# In dir: ~/Projects/Introgressed/genetic_assignment/data
plink --bfile pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo --recode structure --out pruned_merged_domesticated_and_wbdc_318_morex_v3_wPopInfo
```
