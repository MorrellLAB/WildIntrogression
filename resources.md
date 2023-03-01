# Wild barley introgression

- The goals of the project are to identify introgression in barley accessions with a wild phenotype
- There is strong evidence that some of the WBDC accessions have experienced introgression
- Evidence comes from marker-based analysis of similarity to cultivated barley and from phenotypes that aren't fully "wild" such as changes in the degree of shattering.
- Given the many domestication-related genes in barley, chromosome-level genetic assignment could identify regions where domestication related alleles have or have not been introgressed into wild barley.

## Genotyping data

- [Fang et al. 2014][1] report analysis of the [BOPA][2] SNP set in the WBDC collection.
- The full genotyped dataset is available as raw .tsv with Alchemy genotype calls are in the [DRUM][3]
- The WBDC has also been genotyped with GBS as reported in [Sallam et al. 2017][4]
    - The full, filtered dataset is at this path: `/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/Ahmad_GBS_morex_v3`
    - The file with 314 individuals is `final_filtered_treatmissing_biallilic_all_chrs_mafgt0.0016.recode.vcf` but it may be necessary to further filter variants with high heterozygosity or other issues.
        - Additional filtering may be needed to match filtering in the file with only 307 individuals: `WBDC_GBS_noIntro_307_mafgt0.0016.recode.tags.vcf.gz`
- The iSelect 9K dataset from the USDA Core Collection includes most, but not all the BOPA SNPs. A large set of cultivated and landrace lines were genotyped and reported by [Muñoz-Amatriaín et al. 2014][5]
    - The dataset as reported by [Poets et al. 2015][6] is available on [Figshare][7]. There is a GitHub [repository][8] that may be useful for processing this dataset.

## Resequencing data

The fully processed exome capture data set including the introgressed individuals will be in this directory once filtered:

```bash
/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered
```

## Genetic assignment should be useful to determine which individuals are introgressed
- [Admixture][9] has a haploid mode that should be useful for comparing WBDC accessions to domesticated samples from the USDA Core Collection.
- It may be necessary to phase the genotyping data for some analyses that follow, especially genetic assignment of chromosomes or chromosome painting.
    - [SHAPEIT4][10] is a current phasing tool that appears to be relatively easy to use.
- [Plink][11] identity by state analysis can identify nearly identical regions. This may still require Plink 1.9.
- [MOSAIC][12] is current tool for chromosome painting. It requires phased data.

## A list of domestication-related genes is available

- The genes most associated with domestication in barley involve loss of shattering, loss of seed dormancy, and row number. At least one gene has been cloned for each of these phenotypes. We may need to identify Morex_v3 positions.
- The file is [here][13].

[1]: https://doi.org/10.1534/g3.114.010561
[2]: https://doi.org/10.1186/1471-2164-10-582
[3]: https://doi.org/10.13020/D6B59N
[4]: https://doi.org/10.1534/g3.117.300222
[5]: https://doi.org/10.1371/journal.pone.0094688
[6]: https://doi.org/10.1186/s13059-015-0712-3
[7]: https://figshare.com/articles/dataset/Raw_Genotyping_Data_Barley_landraces_are_characterized_by_geographically_heterogeneous_genomic_origins/1468432
[8]: https://github.com/AnaPoets/BarleyLandraces
[9]: https://dalexander.github.io/admixture/download.html
[10]: https://odelaneau.github.io/shapeit4/
[11]: https://www.cog-genomics.org/plink/1.9/distance
[12]: https://maths.ucd.ie/~mst/MOSAIC/
[13]:
