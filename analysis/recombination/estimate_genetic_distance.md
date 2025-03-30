# Estimate genetic distance (cM) for introgressed segments

## The cM size of introgressed segments provides an estimate of their timing

* We currently have estimates of introgression from FLARE that report the physical size of introgressed regions. The estimates are the distance between queried SNPs in the Morex_v3 genome.
* Recombination rate across the barley genome is variable, so to convert physical distances to cM distance.
* Chaochih provided a genetic map of barley BOPA and 9K genetic positions based on 

## Merging genetic and physical map

* Most relevant should be in the DRUM archive

1. Run [merge_genetic_physical_maps.py](https://github.com/MorrellLAB/WildIntrogression/blob/master/analysis/recombination/merge_genetic_physical_maps.py)
This will create "merged_genetic_physical_map.bed"
* Note that this is a BED file with no actual header.

| Chromosome | Start   | End     | Marker   | Genetic Distance (cM) |
|------------|---------|---------|----------|----------------------|
| chr1H      | 50327   | 50328   | 11_20479 | --                   |
| chr1H      | 50466   | 50467   | 12_10420 | 1.3                  |
| chr1H      | 148231  | 148232  | 11_20373 | --                   |
| chr1H      | 157048  | 157049  | 12_30653 | --                   |
| chr1H      | 376847  | 376848  | 11_21354 | 1.38                 |
...

2. Run [interval_genetic_distance.py](https://github.com/MorrellLAB/WildIntrogression/blob/master/analysis/recombination/interval_genetic_distance.py)
This will create "interval_genetic_distances.tsv" with a summary of the genetic and physical size of each introgression interval.

| Chromosome | Start      | End        | Sample  | Genetic_Distance(cM) | Physical_Distance(Mbp) |
|------------|------------|------------|---------|----------------------|------------------------|
| chr1H      | 19822692   | 20591395   | WBDC016 | 0.36                 | 0.768703               |
| chr1H      | 32536998   | 38014817   | WBDC016 | 2.84                 | 5.477819               |
| chr1H      | 124757009  | 288510386  | WBDC016 | 0.00                 | 163.753377             |
| chr1H      | 330397004  | 333666386  | WBDC016 | 0.00                 | 3.269382               |
| chr1H      | 372381382  | 383482448  | WBDC016 | 1.44                 | 11.101066              |
...

3. Run [plot_interval_lengths.py](https://github.com/MorrellLAB/WildIntrogression/blob/master/analysis/recombination/plot_interval_lengths.py)
This will generate seaborn plots summarizing the size of introgression intervals
