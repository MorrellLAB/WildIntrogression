# hap-ibd

Run hap-ibd (https://github.com/browning-lab/hap-ibd) to get a first pass guess at which wild and domesticated individuals share large segments of HBD/IBD segments. This can help us decide which wild samples to include in our query sample set when running FLARE for local ancestry inference.

### Prepare files

Combine all domesticated and wild (including likely introgressed wild) samples into a single VCF.

```bash
module load bcftools/1.10.2
module load htslib/1.9

# In dir: ~/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed
bcftools merge domesticated_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz wbdc_bopa_snps.polymorphic.filt_miss_het.phased.imputed.vcf.gz -O z -o dom_and_wild_with_introgressed_merged.phased.imputed.vcf.gz
tabix -p vcf --csi dom_and_wild_with_introgressed_merged.phased.imputed.vcf.gz

# Exclude SNPs with missing genotypes
bcftools view -e "GT='mis'" dom_and_wild_with_introgressed_merged.phased.imputed.vcf.gz -O z -o dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz
tabix -p vcf --csi dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz
```

### Run hap-ibd

Check average inter-marker distance between SNPs for adjusting hap-ibd parameters.

```bash
# C = prev contig
# P = prev position
# M = max length
# T = total length
# N= num variants
# L = current distance
zcat /panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/vcf/morex_v3/Filtered/phased_and_imputed/dom_and_wild_with_introgressed_merged.phased.imputed.no_missing.vcf.gz | awk 'BEGIN{C="";P=-1;M=0;T=0.0;N=0;} /^#/{next} {if(C==$1) {L=int($2)-P;T+=L;N++;M=(M>L?M:L);}C=$1;P=int($2);}END{if(N>0) printf("avg=%f max=%d\n",T/N,M);}'
avg=2394625.459816 max=125448295
```

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/hap-ibd
sbatch run_hap-ibd-wbdc_nsgc_geno.sh
```

Pull out only pairs with WBDC samples since that's what we're primarily interested in.

```bash
# In dir: ~/Projects/Introgressed/hap_ibd/min2markers
#zgrep "WBDC" dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz > wild_dom_pairs_only.hap-ibd.out.ibd

cd ~/Projects/Introgressed/hap_ibd/minmarkers30_maxgap1000
zgrep "WBDC" dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz | awk '($1~/^WBDC/ && $2!~/^WBDC/) || ($1!~/^WBDC/ && $2~/^WBDC/) { print $0}'
```

### Visualize hap-ibd output

```bash
module load R/4.1.0

cd ~/Projects/Introgressed/hap_ibd/minmarkers50_maxgap1000
~/GitHub/WildIntrogression/analysis/hap-ibd/plot_wbdc_hap-ibd_size_dist.R dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz ~/Projects/Introgressed/hap_ibd/plots minmarkers50_maxgap1000

cd ~/Projects/Introgressed/hap_ibd/minmarkers40_maxgap1000
~/GitHub/WildIntrogression/analysis/hap-ibd/plot_wbdc_hap-ibd_size_dist.R dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz ~/Projects/Introgressed/hap_ibd/plots minmarkers40_maxgap1000

cd ~/Projects/Introgressed/hap_ibd/minmarkers30_maxgap1000
~/GitHub/WildIntrogression/analysis/hap-ibd/plot_wbdc_hap-ibd_size_dist.R dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz ~/Projects/Introgressed/hap_ibd/plots minmarkers30_maxgap1000

cd ~/Projects/Introgressed/hap_ibd/minmarkers20_maxgap1000
~/GitHub/WildIntrogression/analysis/hap-ibd/plot_wbdc_hap-ibd_size_dist.R dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz ~/Projects/Introgressed/hap_ibd/plots minmarkers20_maxgap1000

cd ~/Projects/Introgressed/hap_ibd/minmarkers10_maxgap1000
~/GitHub/WildIntrogression/analysis/hap-ibd/plot_wbdc_hap-ibd_size_dist.R dom_and_wild_with_introgressed_geno.hap-ibd.out.ibd.gz ~/Projects/Introgressed/hap_ibd/plots minmarkers10_maxgap1000
```

Check number of WBDC individuals with IBD segments from domesticated accessions:

```bash
cd ~/Projects/Introgressed/hap_ibd/plots

cut -f 3 minmarkers20_maxgap1000.wbdc-dom_pairs_only.txt | grep "WBDC" | uniq | sort -uV | wc -l
318
cut -f 3 minmarkers30_maxgap1000.wbdc-dom_pairs_only.txt | grep "WBDC" | uniq | sort -uV | wc -l
302
cut -f 3 minmarkers40_maxgap1000.wbdc-dom_pairs_only.txt | grep "WBDC" | uniq | sort -uV | wc -l
51
cut -f 3 minmarkers50_maxgap1000.wbdc-dom_pairs_only.txt | grep "WBDC" | uniq | sort -uV | wc -l
12

cut -f 3 minmarkers50_maxgap1000.wbdc-dom_pairs_only.txt | grep "WBDC" | uniq | sort -uV
WBDC016
WBDC020
WBDC045
WBDC053
WBDC172
WBDC173
WBDC182
WBDC184
WBDC280
WBDC283
WBDC311
WBDC344
```

We'll start with the 40 markers cutoff and use the list of WBDC individuals identified there as the query samples for FLARE.

```bash
cut -f 3 minmarkers40_maxgap1000.wbdc-dom_pairs_only.txt | grep "WBDC" | uniq | sort -uV > potentially_introgressed_wbdc_list.minmarkers40_maxgap1000.txt
```
