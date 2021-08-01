import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab
from studies import targets, sample_timepoints
import re

studies = snakemake.params.studies
N_studies = len(studies)
genes_ID = snakemake.params.genes_ID
genes_symbol = snakemake.params.genes_symbol

READ_UNITS = 1e6 # Use millions of reads
DPI = 300

counter = 0

for gene in genes_ID: 
    fig, ax = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
    fig1, ax1 = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
    #fig2, ax2 = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
    name = genes_symbol[counter]
    for study in studies:
        time_data = []
        time_mod24_data = []
        expression_data = []
        times = sample_timepoints(study)
        times_mod24 = [t % 24 for t in times]
        print(times)
        time_data.extend(times)
        time_mod24_data.extend(times_mod24)
        tpm = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
        for col in tpm.columns:
            expression_data.append(tpm.loc[gene][col])
        ax.scatter(time_data,expression_data, label=study, s=17)
        ax1.scatter(time_mod24_data,expression_data, label=study, s=17)
        expression_data_series=pandas.Series(expression_data)
        mean_by_time = expression_data_series.groupby(time_data).mean()
        ax.plot(mean_by_time.index,mean_by_time)
        mean_by_modtime = expression_data_series.groupby(time_mod24_data).mean()
        ax1.plot(mean_by_modtime.index,mean_by_modtime)

    ax.set_xticks(numpy.arange(0, 72, step=3))
    ax.set_xlabel("Time")
    ax.set_ylabel(name + " Expression TPM")
    ax.set_title(gene)
    fig.legend(fontsize = 'x-small')
    fig.tight_layout()
    fig.savefig(snakemake.output[gene], dpi=DPI)

    ax1.set_xticks(numpy.arange(0, 24, step=3))
    ax1.set_xlabel("Time")
    ax1.set_ylabel(name + " Expression TPM")
    ax1.set_title(gene)
    fig1.legend(fontsize = 'x-small')
    fig1.tight_layout()
    fig1.savefig(snakemake.output[gene+"mod24"], dpi=DPI)
    counter +=1
