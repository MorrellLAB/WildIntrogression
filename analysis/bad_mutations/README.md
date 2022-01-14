# BAD_Mutations

Run [BAD_Mutations](https://github.com/MorrellLAB/BAD_Mutations) to predict deleterious variants.

### Loading dependencies

Load dependencies for BAD_Mutations. Run the lines below before calling on the `BAD_Mutations.py` script for interactive processes.

```bash
module load python3/3.6.3_anaconda5.0.1
export PATH=${PATH}:/panfs/roc/groups/9/morrellp/liux1299/Software/BAD_Mutations
```

#### Step 1: Make config

Load dependencies as mentioned above, then generate the config file.

```bash
# Go into our primary working directory
cd /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations
# Generate config file
BAD_Mutations.py setup \
    -b /panfs/roc/groups/9/morrellp/shared/Projects/Mutant_Barley/results/bad_mutations/Genomes \
    -t "Hvulgare" \
    -e 0.05 \
    -c /panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations/config.txt
```

There should now be a file called `config.txt` in the directory `/panfs/roc/groups/9/morrellp/shared/Projects/Introgressed/bad_mutations`.

#### Step 2: Download CDS files

We can skip this step here since we have previously selected and downloaded 72 Angiosperm genomes to use for this analysis. For reference, those steps are documented here: https://github.com/MorrellLAB/Barley_Mutated/tree/master/02_analysis/bad_mutations

#### Step 3: Generate substitutions files

We will run VeP first before we convert the VeP output to the BAD_Mutations substitutions files format.

```bash
./vep_to_subs.sh
```

#### Step 4: Generate alignments and trees

This step overlaps with the Barley Inversions project, so we do not need to re-run it again here. Please refer to the documentation in that project's Github repository for details.

#### Step 5: Predict substitutions

