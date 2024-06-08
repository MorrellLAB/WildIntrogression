# BLAST known domestication related genes

Downloaded the fasta files of each gene from NCBI by providing a list of GenBank IDs in file `GenBankID_list.txt` (from Table S2) to https://www.ncbi.nlm.nih.gov/sites/batchentrez. Then selecting the option to download all fasta files at once. Previously already downloaded a set of FASTA files, pull updated list from that.

From GenBankID_list.txt 1 gene could not be retrieved: MLOC_69611.1. Retrieved a total of 39 out of 40 genes.

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes
# Pull fasta sequences by Genbank ID match (partial string match mode)
seqkit grep -rf ~/GitHub/WildIntrogression/analysis/known_genes/GenBankID_list.txt domestication_genes_gene_sequence.fasta > ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes_2024-06-01/dom_genes.fasta

cd ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes_2024-06-01
grep ">" dom_genes.fasta | cut -d' ' -f 1 | tr '_' '\t' | cut -f 1 | sed 's/>lcl|//' | sort -V | uniq > temp_genbank_ids_in_fasta.txt
diff -y ~/GitHub/WildIntrogression/analysis/known_genes/GenBankID_list_sorted.txt temp_genbank_ids_in_fasta.txt | grep "<" | cut -f 1 > temp_genbank_ids_missing_fasta_seq.txt
```

Did a BLAST (blastn) search using IPK Galaxy (https://galaxy-web.ipk-gatersleben.de/) because they already have the Morex v3 reference genome available in their database there. Selected extended tabular format output. Then filtered with the script below.

BLAST search all genes together except for:

1. KR813335.1 Hordeum vulgare subsp. spontaneum voucher OUH602, partial sequence (Btr1/Btr2)
2. EF067844.1 Hordeum vulgare vrs1 locus, complete sequence; and Hox1 gene, complete cds

Run these separately (pulling the entire fasta sequence by ID returns errors after 1+ hours when doing a BLAST search).

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes_2024-06-01/blast_out
#filter_blastn_known_genes.sh
#filter_blastn_known_genes_cds.sh
OUT_DIR=~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes_2024-06-01/blast_out/filtered

PERCENT_IDT="95"
E_VALUE="0.0001"
BIT_SCORE="60"

~/GitHub/WildIntrogression/analysis/known_genes/filter_blastn.sh blastn_dom_genes.fasta_vs_morex_v3.tabular ${OUT_DIR} ${PERCENT_IDT} ${E_VALUE} ${BIT_SCORE}

~/GitHub/WildIntrogression/analysis/known_genes/filter_blastn.sh blastn_dom_genes_remaining_3.cds.fasta_vs_morex_v3.tabular ${OUT_DIR} ${PERCENT_IDT} ${E_VALUE} ${BIT_SCORE}

~/GitHub/WildIntrogression/analysis/known_genes/filter_blastn.sh megablast_dom_genes.fasta_vs_morex_v3.tabular ${OUT_DIR} ${PERCENT_IDT} ${E_VALUE} ${BIT_SCORE}
```

Combine bed files and merge domestication genes info with gene positions.

```bash
cd ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes_2024-06-01/blast_out/filtered
cat blastn_dom_genes.fasta_vs_morex_v3.top_hits.bed blastn_dom_genes_remaining_3.cds.fasta_vs_morex_v3.top_hits.bed | sort -k1,1 -k2,2n > dom_genes_combined.gene_seq.bed
```

Check the number of known domestication related genes that overlap with introgressed segments.

```bash
# In dir: ~/Dropbox/Projects/Wild_Introgression/Data/domestication_genes
bedtools intersect -a dom_genes_combined.gene_seq.bed -b ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed | uniq | wc -l
      22

bedtools intersect -wa -a dom_genes_combined.gene_seq.bed -b ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed | uniq > ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/known_dom_genes_overlap_introgressed_regions.bed

# Get 50k SNPs in Table S2 that overlap with introgressed segments
cd ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots_updated_chromPlots_our_pos
bedtools intersect -wa -a 50k_markers_dom_genes-no_disease_resistance_genes.bed -b ~/Dropbox/Projects/Wild_Introgression/Analyses/local_ancestry-flare/plots/wbdc_likely_introgressed_segments-breeding.bed | uniq > known_dom_genes_50k_snps_overlap_introgressed_regions.bed
```

```bash
# Get 50k SNP positions for SNPs in Table S2
grep -wf ~/GitHub/WildIntrogression/analysis/known_genes/table_S2_50k_SNPs_uniq.txt ~/GitHub/morex_reference/morex_v3/50k_9k_BOPA_SNP/50k_idt90_noRescuedSNPs.bed > ~/GitHub/WildIntrogression/analysis/known_genes/table_S2_snps_50k_idt90_noRescuedSNPs.bed
```
