# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants. Some steps are shared with other projects, see here for full documentation: https://github.com/MorrellLAB/Barley_Mutated/tree/master/02_analysis/bad_mutations

### Loading dependencies

Load dependencies for BAD_Mutations. Run the lines below before calling on the `BAD_Mutations.py` script for interactive processes.

```bash
module load python3/3.6.3_anaconda5.0.1
source activate /panfs/roc/groups/9/morrellp/liux1299/.conda/envs/bad_mutations
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/BAD_Mutations
```

#### Step 1: Make config

Config was already generated as part of another project, use the following config:

```bash
/panfs/jay/groups/9/morrellp/gfrascar/bad_mutations_scripts/config.txt
```

#### Step 2: Download CDS files

We can skip this step here since we have previously selected and downloaded 72 Angiosperm genomes to use for this analysis. For reference, those steps are documented here: https://github.com/MorrellLAB/Barley_Mutated/tree/master/02_analysis/bad_mutations

#### Step 3: Generate substitutions files

We ran Annovar and converted the output file to BAD_Mutations .subs files. See analysis directory `Annovar` for methods.

#### Step 4: Generate alignments and trees

This step overlaps with the Barley Inversions project, so we do not need to re-run it again here. Please refer to the documentation in that project's Github repository for details (https://github.com/MorrellLAB/Barley_Inversions/tree/master/01_analyses/BAD_Mutations).

#### Step 5: Predict substitutions

Some preparation steps were already run, see `Step 5: Predict substitutions` section from https://github.com/MorrellLAB/Barley_Mutated/tree/master/02_analysis/bad_mutations for details.

Run BAD_Mutations predict. The `bad_mut_predict-dom_and_wild_snps.sh` script stores filepaths and calls on the main script `bad_mut_predict.sh`. The general command to submit arrays is as below, but in practice we submitted them in batches of ~200 array indices for trackability.

**USEFUL:** For each batch of array indices submitted, there will be some array indices that have timed out or failed and need to be re-run. To generate a list of re-run array indices, run the script `get_re-run_array_indices.sh` (TIMEOUT and FAILED), `get_timeout_array_indices.sh` (TIMEOUT only), or `get_failed_array_indices.sh` (FAILED only).

```bash
# In dir: ~/GitHub/WildIntrogression/analysis/bad_mutations
# Full path to per transcript substitutions directory containing .subs files
#	This output is from the VeP_to_Subs.py supporting script or ANNOVAR_to_subs.py script
SUBS_DIR="/panfs/jay/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations/annovar_to_subs/dom_and_wild_snps"
# Sample name will be used as a prefix for outputs
SAMPLE_NAME="dom_and_wild_snps"
# Full path to output directory
OUT_DIR="/scratch.global/liux1299/bad_mutations/predict_output_${SAMPLE_NAME}"
# Full path to a list of primary transcripts, one per line
PRIMARY_TRANSCRIPTS="/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/phytozome13_download_V3_primary_transcript/hvulgare_primary_transcripts_only.txt"

# Prepare subs list that intersects with primary transcripts list
./intersect_primary_transcripts_and_subs.sh ${SUBS_DIR} ${OUT_DIR} ${PRIMARY_TRANSCRIPTS}

# Run bad mutations predict
# Max number of arrays depends on the number of FASTA_LIST_OF_LISTS
sbatch --array=0-209 bad_mut_predict-dom_and_wild_snps.sh
# Re-submit timeout array indices
get_timeout_array_indices.sh 157562199
# JobID 157673093
sbatch --array=0-66,68-152 bad_mut_predict-dom_and_wild_snps.sh
# JobID 157679606
sbatch --array=153-189 bad_mut_predict-dom_and_wild_snps.sh
# JobID 157687553
sbatch --array=190-193 bad_mut_predict-dom_and_wild_snps.sh

get_timeout_array_indices.sh 157673093
# JobID 157734459
sbatch --array=0-6,9-35,37,39-66,68-92,94,96-115,117-150,152 bad_mut_predict-dom_and_wild_snps.sh

get_timeout_array_indices.sh 157679606
# JobID 157804461
sbatch --array=153-165,167-176,180-189 bad_mut_predict-dom_and_wild_snps.sh
get_timeout_array_indices.sh 157687553
# JobID 157804494
sbatch --array=190-193 bad_mut_predict-dom_and_wild_snps.sh
```

Note: lists 194-209 contain chrUn so no output directories are written. Example, see: `/panfs/jay/groups/9/morrellp/shared/Projects/WBDC_inversions/bad_mutations/results/MSA_output/hvulgare_cds_list-194`
