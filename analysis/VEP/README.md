# VeP

Prepare vcf files if they have not been tabix indexed. VeP only works on bgzipped and tabix indexed VCF files.

Run VeP with script.

```bash
sbatch Vep-snps.sh
sbatch Vep-snps_HC.sh
```

**Note:** When running VeP, don't use `--total_length` flag. Turning on this flag messes up the file format for BAD_Mutations `Vep_to_Subs.py` script.

For some reason, the GFF file isn't working well here.
