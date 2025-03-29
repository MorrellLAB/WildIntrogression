# Population divergence using vcftools

* Generate site wise-diversity and Fst for introgressed and non-introgressed samples
* Files are in [Dropbox](https://www.dropbox.com/scl/fo/12y8eh3zkezyi056n02j5/AA-pX-WTKIu43SUX41sTcYs?rlkey=muclohbydtb3nxrcfvsh7mbt7&dl=0); paths are relative to my user space.

```bash
brew install vcftools

cd /Users/pmorrell/Desktop/WBDC_Introgression
VCF=/Users/pmorrell/Desktop/WBDC_Introgression/wbdc_318_BOPA_morex_v3_poly.vcf.gz
NON_INTR=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/nonintrogressed_samples.txt
INTR=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/introgressed_samples.txt 

vcftools --gzvcf ${VCF} --site-pi --keep ${NON_INTR} --out nonintrogressed

vcftools --gzvcf ${VCF} --site-pi --keep  ${INTR} --out introgressed

vcftools --gzvcf ${VCF} --weir-fst-pop ${NON_INTR} --weir-fst-pop ${INTR}  --out pairwise
```

## `vcftools` is working, now for some analysis

*Create files for cultivated that can be compared to wild introgressed and nonintrogressed

```bash
bcftools view /Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/DRUM_Submission/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz -S ^/Users/pmorrell/Desktop/remove_OWBR.txt | bcftools query -l | sed -e '/WBDC/d' >/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/cultivated_samples.txt 

CULT=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/cultivated_samples.txt
VCF=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/DRUM_Submission/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz

bcftools view /Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/DRUM_Submission/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz -S ^/Users/pmorrell/Desktop/remove_OWBR.txt | bcftools query -l | sed -e '/WBDC/d' >/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/cultivated_samples.txt 

CULT=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/Analyses/introgression_donors/cultivated_samples.txt
VCF=/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Work/Manuscripts/Wild_Introgression/DRUM_Submission/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz
```

* R code for Fst difference between cultivated barley and wild introgressed and nonintrogressed samples

```R
intr = read.table('~/Desktop/WBDC_Introgression/Fst/cult_intr.weir.fst',header=T)
nonintr = read.table('~/Desktop/WBDC_Introgression/Fst/cult_nonintr.weir.fst',header=T)
comp = data.frame(intr[,3],nonintr[,3])
comp_filtered = comp[apply(comp, 1, function(row) all(row > 0)), ]
dim(comp_filtered)
[//]: 1402    2
wilcox.test(comp_filtered[,1],comp_filtered[,2],paired=T)

[//]: 	Wilcoxon signed rank test with continuity correction

[//]: data:  comp_filtered[, 1] and comp_filtered[, 2]
[//]: V = 213768, p-value < 2.2e-16
[//]: alternative hypothesis: true location shift is not equal to 0
```

* Calculate descriptive stats using datamash (sorry R) 

```bash
datamash -H -R 4 --narm mean 3 sstdev 3 median 3 </Desktop/WBDC_Introgression/Fst/cult_intr.weir.fst

datamash -H -R 4 --narm mean 3 sstdev 3 median 3 <Desktop/WBDC_Introgression/Fst/cult_nonintr.weir.fst`
```



## Creating a dxy calculator with Python and GPT4

* vcftools_dxy.py
* **Decided not to use this.** Fst makes more sense here because markers are pre-ascertained so there is not diversity estimate in the dxy comparison.
python ~/Desktop/Sandbox/PMorrell/Utilities/vcftools_dxy.py --pop1 /Users/pmorrell/Desktop/WBDC_Introgression/nonintrogressed.sites.pi --pop2 /Users/pmorrell/Desktop/WBDC_Introgression/introgressed.sites.pi --pairwise /Users/pmorrell/Desktop/WBDC_Introgression/pairwise.weir.fst --output dxy_results.txt
