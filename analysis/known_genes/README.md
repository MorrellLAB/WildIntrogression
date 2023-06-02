# BLAST known domestication related genes

Took a list of `GenBank_ID` from the spreadsheet `` and 

Then, downloaded the fasta files of each gene from NCBI by providing a list of GenBank IDs to https://www.ncbi.nlm.nih.gov/sites/batchentrez. Then selecting the option to download all fasta files at once.

Did a BLAST (blastn) search using IPK Galaxy (https://galaxy-web.ipk-gatersleben.de/) because they already have the Morex v3 reference genome available in their database there. Selected extended tabular format output. Then filtered with the script below.

```bash
filter_blastn_known_genes.sh
filter_blastn_known_genes_cds.sh
```

Combine bed files and merge domestication genes info with gene positions.

```bash
cd ~/Dropbox/Projects/Wild_Introgression/Data/domestication_gene
cat blastn_domestication_genes_gene_sequence.fasta.txt_vs_morex_v3.top_hits.bed vrs1_cds.bed genbank_id_not_in_genes_fasta.cds.sorted.bed | sort -k1,1 -k2,2n > dom_genes_combined.gene_seq.cds.bed
```

Check the number of known domestication related genes that overlap with introgressed segments.

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes
bedtools intersect -a dom_genes_combined.gene_seq.cds.bed -b ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed | uniq | wc -l
      30

bedtools intersect -wa -a dom_genes_combined.gene_seq.cds.bed -b ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed | uniq > ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/known_dom_genes_overlap_introgressed_regions.bed
```
