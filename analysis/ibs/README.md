# Plink IBS/IBD Analysis

We use Plink to perform an IBS analysis (not IBD because we don't provide an outgroup) to detect regions that show introgression.

---

### Summary of Scripts

- `find_z2_threshold.R` script is for data exploration of Plink IBS/IBD output `.genome` files. This script outputs a histogram of the Z2 column and a .txt file with quantiles in 5% intervals of the Z2 column.