# NJ Trees and FLARE Results

Make final figure that combines NJ Trees and FLARE results. Also make publication version of NJ trees for NSGC+WBDC samples and GBS data.

## Methods

### Make list of previously identified wild-introgressed samples

Based on Nice et al. 2016 Table 1 (include Fang et al. 2014 samples identified as likely introgressed).

Number of Fang et al. 2014 wild-introgressed samples identified by FLARE as wild-introgressed.

```bash
wc -l Fang_et_al_2014_wild-introgressed_list.txt
       8 Fang_et_al_2014_wild-introgressed_list.txt

grep -f Fang_et_al_2014_wild-introgressed_list.txt flare_wild-introgressed_table_s3_list.txt | wc -l
       8
```

Nice et al sample list excluding Fang et al identified wild-introgressed samples. Number of these samples identified by FLARE as wild-introgressed.

```bash
grep -vf Fang_et_al_2014_wild-introgressed_list.txt Nice_et_al_2016_25_wild_list.txt > Nice_et_al_2016_exclude_Fang_et_al_2014_wild-introgressed.txt

wc -l Nice_et_al_2016_exclude_Fang_et_al_2014_wild-introgressed.txt
      17 Nice_et_al_2016_exclude_Fang_et_al_2014_wild-introgressed.txt

grep -f Nice_et_al_2016_exclude_Fang_et_al_2014_wild-introgressed.txt flare_wild-introgressed_table_s3_list.txt | wc -l
       6

grep -f Nice_et_al_2016_exclude_Fang_et_al_2014_wild-introgressed.txt flare_wild-introgressed_table_s3_list.txt
WBDC028
WBDC035
WBDC042
WBDC061
WBDC150
WBDC292
```

Half of the samples from this list have slightly longer NJ Tree branches (indicated with "*") in Nice et al. 2016 Figure 1 but are less pronounced as branch lengths for Fang et al. 2014 identified wild-introgressed samples.

WBDC028*
WBDC035
WBDC042*
WBDC061
WBDC150*
WBDC292

### Results Summary

All Fang et al 2014 wild-introgressed samples were also identified as wild-introgressed by FLARE.
