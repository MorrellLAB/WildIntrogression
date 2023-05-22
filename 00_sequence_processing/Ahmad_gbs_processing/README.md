# GBS data filtering

Ahmad's GBS data was filtered as part of the barley inversions project, please see details here: https://github.com/MorrellLAB/Barley_Inversions/tree/master/00_sequence_processing/GBS_Filtering

Li Lei did some filtering (using most of the same criteria as her previous GATK SNP call filtering: https://github.com/MorrellLAB/Barley_Inversions/tree/master/01_analyses/GATK_SNP_call), the file is here:

```bash
/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/Ahmad_GBS_morex_v3/final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.vcf.gz
```

Generate some stats and check proportion het and missingness.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Ahmad_gbs_processing
sbatch check_het.sh
sbatch summarize_missingness.sh
```

Check het/hom_alt ratio for any outliers (generated from bcftools stats output file with some custom calculations in script).

```bash
# In dir: ~/Projects/Introgressed/vcf/Ahmad_GBS_morex_v3/vcf_summary
# het/hom ratios
(head -n 1 final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.het_hom_ratios.txt && tail -n +2 final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.het_hom_ratios.txt | sort -k4,4nr) | head
# Observed heterozygosity
(head -n 1 final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.observed_heterozygosity.txt && tail -n +2 final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.observed_heterozygosity.txt | sort -k5,5nr) | head
```

Check missing.

```bash
sort -k5,5nr final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode_missingness.imiss | head

sort -k6,6nr final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode_missingness.lmiss | head
```

Filter on missingness and heterozygosity.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Ahmad_gbs_processing
sbatch Filter_VCF_Het_and_Miss.sh
```

Re-calculate Ts/Tv.

```bash
# In dir: ~/GitHub/WildIntrogression/00_sequence_processing/Ahmad_gbs_processing
sbatch tstv_snps.sh
```
