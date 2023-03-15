#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=6gb
#SBATCH --tmp=4gb
#SBATCH -t 00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=liux1299@umn.edu
#SBATCH -p small,ram256g,ram1t
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -o pipefail

OUT_DIR="/scratch.global/liux1299/temp_introgression/NSGC_landrace_and_cultivated"
OUT_PREFIX="nsgc_9k"

#------------------
mkdir -p ${OUT_DIR} ${OUT_DIR}/split

cd ${OUT_DIR}/split

# Rawdata_01
wget https://figshare.com/ndownloader/files/2153864
# Rawdata_02
wget https://figshare.com/ndownloader/files/2153865
# Rawdata_03
wget https://figshare.com/ndownloader/files/2153866
# Rawdata_04
wget https://figshare.com/ndownloader/files/2153867
# Rawdata_05
wget https://figshare.com/ndownloader/files/2153868
# Rawdata_06
wget https://figshare.com/ndownloader/files/2153869
# Rawdata_07
wget https://figshare.com/ndownloader/files/2153870
# Rawdata_08
wget https://figshare.com/ndownloader/files/2153874
# Rawdata_09
wget https://figshare.com/ndownloader/files/2153871
# Rawdata_10
wget https://figshare.com/ndownloader/files/2153873
# Rawdata_11
wget https://figshare.com/ndownloader/files/2153872

# Put the pieces together
cat * > ${OUT_DIR}/${OUT_PREFIX}.txt
