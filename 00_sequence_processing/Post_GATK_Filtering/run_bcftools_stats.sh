#!/bin/bash

function generate_stats() {
    local vcf=$1
    local ref_gen=$2
    local out_dir=$3
    # Prep output file prefix
    if [[ "${vcf}" == *".gz"* ]]; then
        prefix=$(basename ${vcf} .vcf.gz)
    else
        prefix=$(basename ${vcf} .vcf)
    fi
    mkdir -p ${out_dir}/plots_${prefix}
    # Generate stats
    bcftools stats -F ${ref_gen} -s - ${vcf} > ${out_dir}/${prefix}.stats
    # Generate plots
    plot-vcfstats -p ${out_dir}/plots_${prefix} -s ${out_dir}/${prefix}.stats
    # Pull out per sample nHets to calculate heterozygosity
    # include header lines
    awk '($1 ~ /PSC/ || $2 ~ /PSC/) { print }' ${out_dir}/${prefix}.stats > ${out_dir}/${prefix}.het_per_sample.txt
    # [4]nRefHom
    # [5]nNonRefHom
    # [6]nHets
    # nHet (RA) / nHom (AA)
    # Calculate het/hom ratio (the value showed in bcftools stats output plot/pdf)
    printf "#[3]sample\t[5]nNonRefHom\t[6]nHets\thet_hom_ratio\n" > ${out_dir}/${prefix}.het_hom_ratios.txt
    grep -v "#" ${out_dir}/${prefix}.het_per_sample.txt | awk '{ OFS="\t" } {$15=$6/$5; print $3,$5,$6,$15}' >> ${out_dir}/${prefix}.het_hom_ratios.txt
    # Calculate observed heterozygosity
    printf "#[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\tobserved_heterozygosity\n" > ${out_dir}/${prefix}.observed_heterozygosity.txt
    grep -v "#" ${out_dir}/${prefix}.het_per_sample.txt | awk '{ OFS="\t" } {$15=$6/($4+$5+$6); print $3,$4,$5,$6,$15}' >> ${out_dir}/${prefix}.observed_heterozygosity.txt
}

export -f generate_stats
