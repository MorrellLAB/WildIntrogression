# Introgression segment length and age estimation

cM/Mb

For Morex_v3 BOPA and 9k SNPs, use Li's script from Environmental Associations study:

- https://github.com/MorrellLAB/Env_Assoc/blob/master/script/recombination/Calculate_cMMb.py
- https://github.com/MorrellLAB/Env_Assoc/blob/master/script/recombination/Smooth_cMMb.R

BOPA and 9k SNP morex_v3 positions taken from: https://github.com/MorrellLAB/morex_reference/tree/master/morex_v3/50k_9k_BOPA_SNP

Prepare VCFs.

```bash
# In dir: ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP
# Concatenate and remove duplicate snps
bcftools concat -a --rm-dups snps bopa_idt90_noRescuedSNPs.vcf.gz 9k_idt90_noRescuedSNPs.vcf.gz | bcftools norm --targets "^chrUn" -d snps -O v -o ~/Dropbox/Projects/Wild_Introgression/Analyses/recombination/bopa_idt90_and_9k_idt90_noRescuedSNPs.vcf

# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/recombination
# Check remaining duplicate SNPs that didn't get removed
grep -v "#" bopa_idt90_and_9k_idt90_noRescuedSNPs.vcf | cut -f 3 | sort -V | uniq -c | awk '$1!="1" { print }'
```

Prepare genetic map files.

```bash
# In dir: ~/GitHub/morex_reference/genetic_maps
awk '{ print $1,$3,$4 }' 9k_iSelect_Genetic_Map.txt | tr ' ' '\t' | tail -n +2 > 9k_iSelect_Genetic_Map_pos.txt
awk '{ print $2,$1,$3 }' BOPA/2013_iSelect_Genetic_Map.txt | tr ' ' '\t' > BOPA/bopa_genetic_map_pos.txt

# Concatenate and remove duplicates
# Check first if duplicates remain
cat BOPA/bopa_genetic_map_pos.txt 9k_iSelect_Genetic_Map_pos.txt | sort -k2,2 -k3,3n | uniq | cut -f 1 | sort -V | uniq -c | awk '$1!="1" { print }'
# Save
cat BOPA/bopa_genetic_map_pos.txt 9k_iSelect_Genetic_Map_pos.txt | sort -k2,2 -k3,3n | uniq > ~/Dropbox/Projects/Wild_Introgression/Analyses/recombination/bopa_and_9k_genetic_pos.txt
```

Run Li's scripts:

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/recombination
# Note: script is for Python v2
python ~/GitHub/Env_Assoc/script/recombination/Calculate_cMMb.py bopa_idt90_and_9k_idt90_noRescuedSNPs.with_duplicates.vcf bopa_and_9k_genetic_pos.txt > bopa_9k_cM_100kb_unsorted.txt
# Sort
head -n 1 bopa_9k_cM_100kb_unsorted.txt > bopa_9k_cM_100kb.txt
tail -n +2 bopa_9k_cM_100kb_unsorted.txt | sort -k1,1 -k2,2n >> bopa_9k_cM_100kb.txt
```

Run a modified version of Li's script: https://github.com/MorrellLAB/Env_Assoc/blob/master/script/recombination/Smooth_cMMb.R
