import json
import re
import pathlib

import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

from studies import targets, sample_timepoints
from styles import color_by_study, shape_by_study, linestyle_by_study

studies = snakemake.params.studies
N_studies = len(studies)
genes_ID = snakemake.params.genes_ID
genes_symbol = snakemake.params.genes_symbol

READ_UNITS = 1e6 # Use millions of reads
DPI = 300


# Load the data
tpm_data = {}
times_data = {}
for study in studies:
    tpm_data[study] = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
    times_data[study] = sample_timepoints(study)

max_time = max(t for times in times_data.values() for t in times)


# Plot genes two different ways:
# - By collapsing periods onto each other by mod 24
# - without collapsing periods onto each other
for name, gene in zip(genes_symbol, genes_ID): 
    fig, ax = pylab.subplots(figsize=(2+max_time/6,4))
    fig1, ax1 = pylab.subplots(figsize=(2+max_time/6,4))
    for study in studies:
        time_data = []
        time_mod24_data = []
        expression_data = []
        times = times_data[study]
        times_mod24 = [t % 24 for t in times]
        time_data.extend(times)
        time_mod24_data.extend(times_mod24)
        tpm = tpm_data[study]
        for col in tpm.columns:
            expression_data.append(tpm.loc[gene][col])
        ax.scatter(time_data,expression_data, label=study, s=17, marker=shape_by_study[study], color=color_by_study[study])
        ax1.scatter(time_mod24_data,expression_data, label=study, s=17, marker=shape_by_study[study], color=color_by_study[study])
        expression_data_series=pandas.Series(expression_data)
        mean_by_time = expression_data_series.groupby(time_data).mean()
        ax.plot(mean_by_time.index,mean_by_time, color=color_by_study[study], linestyle=linestyle_by_study[study])
        mean_by_modtime = expression_data_series.groupby(time_mod24_data).mean()
        ax1.plot(mean_by_modtime.index,mean_by_modtime, marker=shape_by_study[study], color=color_by_study[study])

    ax.set_xticks(numpy.arange(0, 72, step=6))
    ax.set_xticks(numpy.arange(0, 72, step=3), minor=True)
    ax.set_xlabel("Time")
    ax.set_ylabel(name + " Expression TPM")
    ax.set_title(f"{gene} | {name}")
    ax.margins(x=0.01)
    fig.legend(fontsize = 'x-small', ncol=2)
    fig.tight_layout()
    fig.savefig(snakemake.output[gene], dpi=DPI)

    ax1.set_xticks(numpy.arange(0, 24, step=6))
    ax1.set_xticks(numpy.arange(0, 24, step=3), minor=True)
    ax1.set_xlabel("Time")
    ax1.set_ylabel(name + " Expression TPM")
    ax1.set_title(f"{gene} | {name}")
    ax1.margins(x=0.01)
    fig1.legend(fontsize = 'x-small',ncol=2)
    fig1.tight_layout()
    fig1.savefig(snakemake.output[gene+"mod24"], dpi=DPI)
