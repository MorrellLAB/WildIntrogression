# Draw trees to inform sample selection

Draw trees to decide on an intial set of WBDC samples to include in query panel for looking at potential introgression from domesticated samples.

### Convert VCF to PHYLIP format

```bash
# In dir: ~/Software/vcf2phylip
module load python3/3.8.3_anaconda2020.07_mamba

OUT_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees"
```

```bash
# WBDC only BOPA/9k genotypes
./vcf2phylip.py --input /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/wbdc_bopa_snps.polymorphic.filt_miss_het.vcf.gz --output-folder ${OUT_DIR}
```

```bash
# NSGC and WBDC BOPA/9k genotypes
./vcf2phylip.py --input /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz --output-folder ${OUT_DIR}
# The OWBDR accession names are >10 characters long which results in PHYLIP shortening the names so they are no longer unique
# which causes a "Error: duplicate taxon OWBDR60240" error message in figtree
# We'll give these accessions a different name that we can use a lookup table later
grep OWB nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phy | awk '{ print $1 }' > nsgc_names_too_long_for_phylip.txt
# Create a lookup of full accession name and short name
seq 1 29 > tmp_1-29.txt
sed -i 's,^,OWBDR,' tmp_1-29.txt
paste nsgc_names_too_long_for_phylip.txt tmp_1-29.txt | sed 's,\t, ,' > nsgc_short_names_lookup_for_OWBDR.txt

# Rename samples in VCF
module load bcftools/1.10.2
bcftools reheader --samples nsgc_short_names_lookup_for_OWBDR.txt /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.vcf.gz > nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.short_names.vcf.gz
# module load python3/3.8.3_anaconda2020.07_mamba
# ~/GitHub/WildIntrogression/analysis/trees/replace_long_names_phylip_outfile.py nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phy nsgc_short_names_lookup_for_OWBDR.txt > nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.renamed_samp.phy

# Re-run vcf2phylip
./vcf2phylip.py --input /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.short_names.vcf.gz --output-folder ${OUT_DIR}
```

```bash
# WBDC GBS data
./vcf2phylip.py --input /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/Ahmad_GBS_morex_v3/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.vcf.gz --output-folder ${OUT_DIR}
```

```bash
# Domesticated and WBDC exome capture
# This one takes longer and needs to be submitted as a Slurm job script
sbatch run_vcf2phylip-exome_cap.sh
```

### Generate distance matrices with PHYLIP (PHYLogeny Inference Package)

Use [PHYLIP](https://evolution.genetics.washington.edu/phylip) to generate distance matrices for drawing trees.

```bash
# Add to path
export PATH=${PATH}:/panfs/roc/msisoft/phylip/3.69/bin/

# Go into working directory
cd ~/Projects/Introgressed/trees
# Invoke program
dnadist
# Prompt for input file, paste the following:
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/wbdc_bopa_snps.polymorphic.filt_miss_het.min4.phy
# Menu options will pop up
# Accept defaults
Y
# Rename output file
mv outfile wbdc_bopa_snps.polymorphic.filt_miss_het.min4.phylip.outfile

# Run neighbor
neighbor
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/wbdc_bopa_snps.polymorphic.filt_miss_het.min4.phylip.outfile
Y
# Rename output files
mv outtree wbdc_bopa_snps.polymorphic.filt_miss_het.min4.phylip.outtree
mv outfile wbdc_bopa_snps.polymorphic.filt_miss_het.min4.phylip.outtree.outfile
```

Visualize `.outtree` file with figtree: https://github.com/rambaut/figtree

Other datasets are larger and get terminated when run on MSI's login nodes, so we'll submit a Slurm interactive job to run the rest.

```bash
srun -N 1 --ntasks-per-node=8 --mem=12gb --tmp=8gb -t 2:00:00 -p interactive --pty bash

# Once job has started, run the following
cd ~/Projects/Introgressed/trees
export PATH=${PATH}:/panfs/roc/msisoft/phylip/3.69/bin/

#----------
# GBS data
# Invoke program
dnadist
# Prompt for input file, paste the following:
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.min4.phy
# Menu options will pop up
# Accept defaults
Y
# Rename output file
mv outfile WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.min4.phylip.outfile

# Run neighbor
neighbor
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.min4.phylip.outfile
Y
mv outtree WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.min4.phylip.outtree
mv outfile WBDC_GBS_snps_morex_v3.biallelic.mafgt0.0016.filt_miss_het.min4.phylip.outtree.outfile

#----------
# Genotype NSGC and WBDC
# Invoke program
dnadist
# Prompt for input file, paste the following:
#/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phy
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.short_names.min4.phy
# Menu options will pop up
# Accept defaults
Y
# Rename output file
mv outfile nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phylip.outfile

# Run neighbor
neighbor
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phylip.outfile
Y
mv outtree nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phylip.outtree
mv outfile nsgc_wbdc_bopa_9k_morex_v3.poly.filt_miss_het.min4.phylip.outtree.outfile

#----------
# Exome capture data
# Invoke program
dnadist
# Prompt for input file, paste the following:
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/dom_and_wild_snps_biallelic.callable.cap50x.final.min4.phy
# Menu options will pop up
# Accept defaults
Y
# Got the error message:
ERROR: bad base: E at site     1 of species 545

# Try removing first sample
grep -v "100-HS8" dom_and_wild_snps_biallelic.callable.cap50x.final.min4.phy > dom_and_wild_snps_biallelic.callable.cap50x.final.min4.excluded_samp.phy

# We'll try running dnadist again
dnadist
/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/trees/dom_and_wild_snps_biallelic.callable.cap50x.final.min4.excluded_samp.phy
# Still get the same error
ERROR: bad base: E at site     1 of species 544

# Likely an issue with the original file conversion
```


