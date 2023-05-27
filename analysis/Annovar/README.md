# ANNOVAR

Run [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) to annotate intergenic, exonic, intronic, synonymous, and nonsynonymous SNPs.

**Purpose:** To prepare the nonsynonymous SNPs for predicting deleterious SNPs using BAD_Mutations and other downstream tools.

### Files required

- VCF file: `/panfs/jay/groups/9/morrellp/shared/Projects/Mutant_Barley/de_novo_vcfs/mut8_and_3mut10xGenomics.SNPs.private.vcf.gz`
- Reference genome file: `/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta`
- GFF3 file(s):
    - All: `/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.parts.gff3`

### ANNOVAR Steps

We will mostly follow the guide here: https://annovar.openbioinformatics.org/en/latest/user-guide/gene/. We'll run it for both high confidence gene models and all gene models.

Load dependencies for ANNOVAR.

```bash
module load perl/5.26.1
# Export paths to directories containing annovar scripts so scripts are calleble from anywhere without specifying the path
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar
export PATH=${PATH}:/panfs/jay/groups/9/morrellp/shared/Software/annovar_conversion_tools
```

---

### Run Annovar on gene models

Filepath

```bash
#GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/PhytozomeV13_HvulgareMorex_V3/annotation/HvulgareMorex_702_V3.gene.parts.gff3"
GFF="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.sorted.parts.gff3.gz"

OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar"

REF_FASTA="/panfs/jay/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex_v3/Barley_MorexV3_pseudomolecules_parts.fasta"

VCF="/panfs/jay/groups/9/morrellp/shared/Datasets/Alignments/introgression_project/all_dom_and_wild/Filtered/dom_and_wild_snps_biallelic.callable.cap50x.final.vcf.gz"
```

Generate the genePred file from the GFF3 file.

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/Annovar
./gff3_to_gene_pred.sh ${GFF} ${OUT_DIR} hv_morex_v3_hc
```

Build the database: generate a transcript FASTA file.

```bash
./build_db.sh ${OUT_DIR}/hv_morex_v3_hc_refGene.txt ${REF_FASTA} ${OUT_DIR} hv_morex_v3_hc
```

Prepare VCF file for Annovar. Convert VCF to Annovar's input format using Annovar's `convert2annovar.pl` script.

```bash
#./vcf_to_annovar_input.sh ${VCF} ${OUT_DIR}
sbatch vcf_to_annovar_input.sh
```

Annotation with Annovar using the script `annotate_variation.pl`.

```bash
#./annotate_with_annovar.sh "/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/Annovar/dom_and_wild_snps_biallelic.callable.cap50x.final_annovar_input.txt" ${OUT_DIR} "hv_morex_v3_hc"
sbatch annotate_with_annovar.sh
```

Description of nomenclature in column 3 of `.exonic_variant_function`: https://varnomen.hgvs.org/bg-material/simple/

Quick exploration of Annovar output.

```bash
# In dir: ~/Projects/Introgressed/Annovar
cut -f 2 dom_and_wild_snps_biallelic.callable.cap50x.final_annovar_input.txt.exonic_variant_function | sort -uV
nonsynonymous SNV
stopgain
stoploss
synonymous SNV
unknown

liux1299@ln0003:~/Projects/Introgressed/Annovar $ cut -f 1 dom_and_wild_snps_biallelic.callable.cap50x.final_annovar_input.txt.variant_function | sort -uV
UTR3
UTR5
UTR5;UTR3
downstream
exonic
exonic;splicing
intergenic
intronic
splicing
upstream
upstream;downstream
```

---

### Prepare Annovar output for BAD_Mutations

Convert Annovar output file to BAD_Mutations .subs format.

```bash
sbatch annovar_to_subs.sh
```

Check total number of .subs files.

```bash
# In dir: ~/Projects/Introgressed/bad_mutations/annovar_to_subs/dom_and_wild_snps
ls *.subs | wc -l
21700
```
