# Meta-analysis of diurnal transcriptomics reveals strong patterns of concordance and discordance in mouse liver
Comparison of 'control' circadian/diurnal timeseries data across multiple studies

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
