import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab
import statsmodels.api as sm
import util
import scipy.stats
from studies import sample_timepoints

studies = snakemake.params.studies
DPI = 300

robustness_score = pandas.read_csv(snakemake.input.robustness, sep="\t", index_col=0)["0"]
data = {}
for study, tpmfile in zip(studies, snakemake.input.tpm):
    if study == "Weger18":
        continue
    tpm = pandas.read_csv(tpmfile, sep="\t", index_col=0)
    data[study] = tpm

data_df = pandas.concat(data.values(), axis=1)

med_tpm = data_df.median(axis=1)

fig, ax = pylab.subplots(figsize=(12,12))
ax.scatter(med_tpm, robustness_score+numpy.random.normal(size=len(robustness_score))*0.2)
ax.set_xscale("log")
ax.set_xlabel("Median TPM")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.expression_level, dpi=DPI)

#amp, period, phase from jtk, gc content

data = {}
for study, ampfile in zip(studies, snakemake.input.jtk):
    if study == "Weger18":
        continue
    amp = pandas.read_csv(ampfile, sep="\t", index_col=0)["AMP"]
    data[study] = amp

data_df = pandas.concat(data.values(), axis=1)

med_amp = data_df.median(axis=1)

fig, ax = pylab.subplots(figsize=(12,12))
ax.scatter(med_amp, robustness_score+numpy.random.normal(size=len(robustness_score))*0.2)
ax.set_xscale("log")
ax.set_xlabel("Median Amplitude")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.amplitude, dpi=DPI)

data = {}
for study, perfile in zip(studies, snakemake.input.jtk):
    if study == "Weger18":
        continue
    per = pandas.read_csv(perfile, sep="\t", index_col=0)["PER"]
    data[study] = per

data_df = pandas.concat(data.values(), axis=1)

med_per = data_df.median(axis=1)
#mean_per = scipy.stats.circmean(data_df, low=0, high=24, axis=1)

fig, ax = pylab.subplots(figsize=(12,12))
ax.scatter(med_per+numpy.random.normal(size=len(robustness_score))*0.2, robustness_score+numpy.random.normal(size=len(robustness_score))*0.2)
#ax.set_xscale("log")
ax.set_xlim(19,29)
ax.set_xlabel("Median Period")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.period, dpi=DPI)

data = {}
for study, phasefile in zip(studies, snakemake.input.jtk):
    if study == "Weger18":
        continue
    jtk = pandas.read_csv(phasefile, sep="\t", index_col=0)
    p=jtk["ADJ.P"]
    phase=jtk["LAG"]
    phase[p>0.05]=float("nan")
    data[study] = phase

data_df = pandas.concat(data.values(), axis=1)

#med_phase = data_df.median(axis=1)
mean_phase = scipy.stats.circmean(data_df, low=0, high=24, axis=1, nan_policy="omit")

fig, ax = pylab.subplots(figsize=(12,12))
ax.scatter(mean_phase, robustness_score+numpy.random.normal(size=len(robustness_score))*0.2)
#ax.set_xscale("log")
ax.set_xlabel("Mean Phase")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.phase, dpi=DPI)
