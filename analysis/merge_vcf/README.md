# Merging VCF files

The `merge_vcf.py` script takes in two vcf files and a contig file containing `##contig` positions (see below for reason). This script outputs an unsorted and NOT de-duplicated merged VCF file, so two steps are needed to get a sorted uniq VCF file:

1. `./merge_vcf.py toy_domesticated.vcf toy_318_wbdc_BOPA.vcf contig_header.vcf > test_merged.vcf`
2. `(head -n 18 test_merged.vcf && tail -n +19 test_merged.vcf | sort -uV -k1,2) > test_merged_sorted.vcf`

Where `head -n 18` refers to the number of meta-info (`##`) and header (`#CHROM`) lines. We then do an alphanumeric unique sort on columns 1 and 2 (#CHROM and POS columns respectively) before saving to an output file.

---

#### Reason we need to provide a file with contig lengths

Given the contig positions output from Plink is inaccurate (because Plink calculates contig length based on last SNP), the `merge_vcf.py` script replaces the contig lines with WBDC exome capture contig lengths. As a test, I will use [`bcftools`](https://samtools.github.io/bcftools/bcftools.html) on toy data to work out meta-info merge process. This ensures my merged VCF complies with the standard VCF file format.

NOTE: the reason we are **NOT** using `bcftools` to merge full length VCF files is because there is a `.bai` limitation with maximum index number at 2^29 - 1.
