# Verify SNPs in Alchemy genotype data

We'll mostly follow the tutorial here for fixing reference strand issues in the BOPA genotype data: https://github.com/MorrellLAB/Barley_Inversions/blob/master/01_analyses/SNP_valiadation/tutorial_alchemy2vcf.md

---

We'll skip the original "Step 1: Get Plink 1.9 ped and map files" and "Step 2: Run Plink 1.9 to recode map and ped file to vcf file" since we already have a VCF file to test.

```bash
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf
```

### Prepare VCF file and other required files

Not all BOPA SNPs are included in the 9k SNPs VCF. So we'll concatenate the two and keep unique SNPs so the following check are easier. The NSGC and WBDC datasets were genotyped on BOPA and 9K.

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP
module load htslib/1.9
# Prepare BOPA and 9k VCFs
bgzip -c bopa_idt95_noRescuedSNPs.vcf > ~/scratch/bopa_idt95_noRescuedSNPs.vcf.gz
bgzip -c 9k_idt95_noRescuedSNPs.vcf > ~/scratch/9k_idt95_noRescuedSNPs.vcf.gz
tabix -p vcf --csi ~/scratch/bopa_idt95_noRescuedSNPs.vcf.gz
tabix -p vcf --csi ~/scratch/9k_idt95_noRescuedSNPs.vcf.gz
# Concatenate BOPA and 9k VCFs
bcftools concat --rm-dups snps --allow-overlaps ~/scratch/bopa_idt95_noRescuedSNPs.vcf.gz ~/scratch/9k_idt95_noRescuedSNPs.vcf.gz | bcftools sort -O v -o ~/scratch/temp_bopa_9k_idt95_noRescuedSNPs.vcf
# Check for duplicates
tail -n +20 ~/scratch/temp_bopa_9k_idt95_noRescuedSNPs.vcf | sort -k1,1 -k2,2n | cut -f 3 | sort -V | uniq -c | sort -k1,1r | head
# Sort unique
(head -n 19 ~/scratch/temp_bopa_9k_idt95_noRescuedSNPs.vcf && tail -n +20 ~/scratch/temp_bopa_9k_idt95_noRescuedSNPs.vcf | sort -k1,1 -k2,2n) | uniq > ~/scratch/bopa_9k_idt95_noRescuedSNPs.vcf
```

Prepare output directory.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3
mkdir fix_strand
cd fix_strand

cp /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/merged_domesticated_and_wbdc_318_morex_v3.vcf test.vcf
```

Load dependencies.

```bash
module load vcftools_ML/0.1.16
module load bcftools/1.10.2
```

The `test.vcf` file includes the following categories of solutions:

1. Force reference
2. Flipped strand
3. First flipped and then force reference

We'll use vcftools v0.1.16 to figure out what needs to be fixed.

First, sort VCF file.

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/fix_strand
vcf-sort test.vcf > sorted_test.vcf
```

Then, run vcftools to figure out the concordant and discordant SNPs:

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/fix_strand
vcftools --vcf sorted_test.vcf --diff ~/scratch/bopa_9k_idt95_noRescuedSNPs.vcf --diff-site --out diff_test_BOPA_9k
```

---

### Step 3: Find concordant and discordant SNPs

- Concordant SNPs with BOPA - "B"
- Only called in BOPA file - "2"
- Discordant SNPs with BOPA - "O"

```bash
# In dir: ~/Projects/Introgressed/vcf/morex_v3/fix_strand
grep "B" diff_test_BOPA_9k.diff.sites_in_files | wc -l
    2963

awk '$4=="2"{print $0}' diff_test_BOPA_9k.diff.sites_in_files | wc -l
    4792

awk '$4=="O"{print $0}' diff_test_BOPA_9k.diff.sites_in_files | wc -l
    0

#awk '$4=="O"{print $0}' diff_test_BOPA_9k.diff.sites_in_files > diff_test_BOPA_9k.discordant.txt
```

Here, we don't have any discordant BOPA/9k SNPs, so we can proceed with the file as is and don't need to run remaining steps.
