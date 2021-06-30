#!/usr/bin/env sh
module load /project/itmatlab/sharedmodules/use.shared
module load sratoolkit-2.10.0
module load salmon-v1.4.0
source venv/bin/activate

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --profile lsf -c 50
