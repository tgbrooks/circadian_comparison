# Meta-analysis of diurnal transcriptomics in mouse liver reveals low repeatability of rhythm analyses
Comparison of 'control' circadian/diurnal timeseries data across multiple studies. This is the pipeline for the paper:

Brooks TG, Manjrekar A, Mrčela A, Grant GR. [Meta-analysis of Diurnal Transcriptomics in Mouse Liver Reveals Low Repeatability of Rhythm Analyses](https://journals.sagepub.com/doi/10.1177/07487304231179600). Journal of Biological Rhythms. 2023.

This project aggregates RNA-seq data from dozens of circadian timeseries transcriptomic studies in the GEO repository.
Datasets were selected to contain mice from 'control' conditions.
While many studies contain disparate conditions, they often contain a timeseries under control conditions.
Aggregating these controls across many studies gives a large number of comparable datasets performed at different time and different institutions.

The data is analyzed for consistency and variation in the rhythmic behavior of the transcripts across these studies.
Any such variation between studies captures technical and biological variance that cannot be observed in any individual study.
This addresses the question of how repeatable can we expect circadian experiments to be: will another lab performing the same experiment get the same list of rhythmic genes?
This further informs comparisons across different conditions, allowing us to compare the variation between conditions to the variations expected just from repeating the same conditions. 

This work began as a [APSA Virtual Summer Research Program](https://www.physicianscientists.org/page/VSRP-2021) in 2021 by Aditi Mandjrekar and mentor Thomas Brooks with Dr. Gregory Grant at the Institute for Translational Medicine and Therapeutics at UPenn.

This repository contains a Snakemake work-flow to download, process, and analyze the RNA-seq data from these datasets.
It is run from `run_pipeline.sh` which also loads necessary software, Snakemake cluster configuration and parameters like job limits.
To run on your cluster, you will likely need to modify this to load the required software packages.

Studies to be processed are listed in `studies.py` under the `targets` dictionary and `studies` list.

## Data Availability

For most users, we recommend the published results of this dataset [available on Zenodo](https://zenodo.org/record/7760579).
This includes tables of quantified values (TPM, read counts) as well as results of running JTK, compareRhythms, booteJTK, compareRhythms and shape invariant models (SIM) across the studies.
This saves considerable effort versus attempting to re-run the pipeline.

## Running the pipeline

The pipeline requires Python 3.9. First, clone this directory and `cd` into it. The python dependencies can be installed by the following.

``` bash
python -m venv venv
source venv/bin/activate # For unix systems, command differs on Windows
pip install -r requirements.txt
```

Additional non-Python dependencies are required to be available on the PATH:

```
sratoolkit-2.11.0
salmon-v1.4.0
R/3.6.3
edirect-15.3
apptainer 1.1.3
```
Next, `R` needs the following libraries installed:

```
tidyverse, MetaCycle, assist, compareRhythms
```

Finally, running BooteJTK is done via `apptainer` to run an image with the appropriate dependencies installed.
First, use `apptainer build` to build the dockerfile `bootejtk/Dockerfile` and then update the directory in `Snakemake` file to point to this (default location is `~/.apptainer/images/ejtk_bootejtk.sif`).

## Snakemake Pipeline

The below image gives the rule dependency graph of the pipeline. Click to enlarge.

![rule dependency graph](https://raw.githubusercontent.com/tgbrooks/circadian_comparison/main/rulegraph_graphviz.svg)

In the Snakefile, some outputs are commented out in the `all` rule (namely, those involving the SIM fits, which involves several hundred hours of compute time).

## Outputs

The pipeline outputs study-level intermediate and processed data to `data/$STUDY_NAME` directories.
Results that aggregate across studies are output to `results/Liver/`.
