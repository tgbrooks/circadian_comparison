#!/usr/bin/env sh
module load /project/itmatlab/sharedmodules/use.shared
module load sratoolkit-2.11.0
module load salmon-v1.4.0
module load R/3.6.3
module load edirect-15.3
source venv/bin/activate


### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --profile lsf -j 50 -c 300
