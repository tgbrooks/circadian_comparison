import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

studies = snakemake.params.studies
N_studies = len(studies)

DPI = 300
Q_CUTOFF = 0.05 # Q-value cutoff to call significantly rhythmic

def num_below(values, cutoffs):
    '''
    Returns the number of entries in `values` that are
    below each cutoff in `cutoffs`
    '''
    values = numpy.sort(values)
    return numpy.searchsorted(values, cutoffs, side="right")

# Breakdowns of the number of q-values less than (or equal to)
# a given cutoff
breakpoints = numpy.linspace(0,1,501) # q-values to compute for
breakdowns = {}
periods = {}
phases = {}
amplitudes = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    breakdowns[study] = num_below(jtk['BH.Q'], breakpoints)
    significant = jtk[jtk['BH.Q'] < Q_CUTOFF]
    periods[study] = significant['PER']
    phases[study] = significant['LAG']
    amplitudes[study] = significant['AMP']

#Q-value breakdowns
fig, ax = pylab.subplots(figsize=(4,4))
for study, bds in breakdowns.items():
    ax.plot(breakpoints, bds, label=study)
ax.set_xlabel("Q-value cutoff")
ax.set_ylabel("Number of Genes")
fig.legend()
fig.tight_layout()
fig.savefig(snakemake.output.breakdowns, dpi=DPI)

# Histogram of periods
fig, axes = pylab.subplots(figsize=(6,0.7+0.7*N_studies), nrows=N_studies, sharex=True)
bins = numpy.linspace(19.5,28.5,10)
for ax, (study, period) in zip(axes.flatten(), periods.items()):
    ax.hist(period, bins=bins)
    ax.set_ylabel(study, rotation=0, horizontalalignment="right")
axes[-1].set_xlabel("Period (hrs)")
fig.tight_layout()
fig.savefig(snakemake.output.periods, dpi=DPI)

# Histogram of phases 
fig, axes = pylab.subplots(figsize=(6,0.7+0.7*N_studies), nrows=N_studies, sharex=True)
bins = numpy.linspace(0,24,25)
for ax, (study, phase) in zip(axes.flatten(), phases.items()):
    ax.hist(phase, bins=bins)
    ax.set_ylabel(study, rotation=0, horizontalalignment="right")
axes[-1].set_xlabel("Phase (hrs)")
fig.tight_layout()
fig.savefig(snakemake.output.phases, dpi=DPI)

# Histogram of amplitudes 
fig, axes = pylab.subplots(figsize=(6,0.7+0.7*N_studies), nrows=N_studies, sharex=True)
bins = numpy.logspace(-1,4, 51)
for ax, (study, amplitude) in zip(axes.flatten(), amplitudes.items()):
    ax.hist(amplitude, bins=bins)
    ax.set_xscale("log")
    ax.set_ylabel(study, rotation=0, horizontalalignment="right")
axes[-1].set_xlabel("Amplitude (TPM)")
fig.tight_layout()
fig.savefig(snakemake.output.amplitudes, dpi=DPI)
