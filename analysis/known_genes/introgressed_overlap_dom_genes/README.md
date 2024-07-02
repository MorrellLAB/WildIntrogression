# Plot Figure S3

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/introgressed_overlap_dom_genes
./prep_bed_files.R

bedtools intersect -wa -a ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed -b dom_genes_50k_snps_noDiseaseRes.bed | sort -k4,4 -k1,1 -k2,2n | uniq > introgressed_segments_overlap_dom_genes_50k_snps_noDiseaseRes.bed

bedtools intersect -wa -a ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed -b dom_genes_noDiseaseRes.bed | sort -k4,4 -k1,1 -k2,2n | uniq > introgressed_segments_overlap_dom_genes_noDiseaseRes.bed

# Get number of unique dom genes that overlap introgressed segment
bedtools intersect -wa -a dom_genes_noDiseaseRes.bed -b ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed | sort -k4,4 -k1,1 -k2,2n | uniq | wc -l
      20
```
