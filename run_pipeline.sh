#!/usr/bin/env sh
module load /project/itmatlab/sharedmodules/use.shared
module load sratoolkit-2.10.0
module load salmon-v1.4.0
source venv/bin/activate


### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --profile lsf -c 50
