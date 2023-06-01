# ABBA-BABA / fd admixture estimation

We'll use the fd statistic described in Martin et al. 2015 MBE that is implemented here: https://github.com/simonhmartin/genomics_general

We will use the same tool used in Wang et al (the citrus introgression paper):
https://github.com/simonhmartin/genomics_general#abba-baba-statistics-in-sliding-windows
- Here's some tutorials:
    - https://github.com/simonhmartin/tutorials/blob/master/ABBA_BABA_whole_genome/README.md
    - https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_windows
    - https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing

We'll have to use large window sizes (chromosome sized) as input to the script because windows that are too small don't work well for the ABBA BABA test because small windows squeeze the data in terms of 
information content.

---
## NSGC and WBDC genotype data

Use the phased and imputed genotype data.

Prepare pops.txt file with 2 columns: 1) sample names and 2) population name

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/ABBA-BABA
module load python3/3.8.3_anaconda2020.07_mamba
#./generate_pops_file.py Sample_Info_Complete.csv vcf_sample_names-dom-wild-Hmurinum_morex_v3.txt > dom-wild-Hmurinum.pops.txt
```

Try only including WBDC individuals with > 10% proportion SNPs assigned to domesticated group.

```bash
cd ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots
grep "breeding" wild_introgressed_prop.txt | sort -k4,4 | awk '$4 >= 0.1 { print }' > wild_introgressed_prop.gte10percent.txt
```

