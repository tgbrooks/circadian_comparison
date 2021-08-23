import json
import math
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

import util

studies = snakemake.params.studies
N_studies = len(studies)

DPI = 300
N_COLUMNS = 3
# Q-value cutoff to call significantly rhythmic
STRICT_Q_CUTOFF = 0.01
LOOSE_Q_CUTOFF = 0.10
Q_CUTOFFS = {
    "loose": LOOSE_Q_CUTOFF,
    "strict": STRICT_Q_CUTOFF,
}
color_by_strictness = {
    "loose": "blue",
    "strict": "darkblue",
}

def num_below(values, cutoffs):
    '''
    Returns the number of entries in `values` that are
    below each cutoff in `cutoffs`
    '''
    values = numpy.sort(values)
    return numpy.searchsorted(values, cutoffs, side="right")

# Breakdowns of the number of q-values less than (or equal to)
# a given cutoff
breakpoints = numpy.logspace(-4,0,501) # q-values to compute for
breakdowns = {}
periods = {strictness:{} for strictness in Q_CUTOFFS.keys()}
phases = {strictness:{} for strictness in Q_CUTOFFS.keys()}
amplitudes = {strictness:{} for strictness in Q_CUTOFFS.keys()}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    breakdowns[study] = num_below(jtk['qvalue'], breakpoints)

    for strictness, cutoff in Q_CUTOFFS.items():
        significant = jtk[jtk['qvalue'] < cutoff]
        periods[strictness][study] = significant['PER']
        phases[strictness][study] = significant['LAG']
        amplitudes[strictness][study] = significant['AMP']


#Q-value breakdowns
fig, ax = pylab.subplots(figsize=(4,4))
for study, bds in breakdowns.items():
    ax.plot(breakpoints, bds, label=study)
    ax.set_xscale("log")
ax.set_xlabel("Q-value cutoff")
ax.set_ylabel("Number of Genes")
fig.legend(fontsize = 'x-small')
fig.tight_layout()
fig.savefig(snakemake.output.breakdowns, dpi=DPI)

# Exclude any studies from plotting if they have little significant at the looser cutoff
include_studies = [study for study in studies if len(phases['loose'][study]) > 10]
N_studies = len(include_studies)
N_ROWS = math.ceil(N_studies / N_COLUMNS)
remove_unplotted = slice(0,0) if N_studies == N_ROWS * N_COLUMNS else slice(-(N_ROWS * N_COLUMNS) - N_studies, None)

# Histogram of periods
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.linspace(19.5,28.5,10)
for strictness, periods_ in periods.items():
    for ax, study in zip(axes.flatten(), include_studies):
        period = periods_[study]
        ax.hist(period, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study, rotation=0, horizontalalignment="right")
[ax.set_xlabel("Period (hrs)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(snakemake.output.periods, dpi=DPI)

# Histogram of phases
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.linspace(0,24,25)
for strictness, phases_ in phases.items():
    for ax, study in zip(axes.flatten(), include_studies):
        phase = phases_[study]
        ax.hist(phase, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study, rotation=0, horizontalalignment="right")
        ax.set_xlim(0,24)
        ax.set_xticks(numpy.arange(0,30, 6))
[ax.set_xlabel("Phase (hrs)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(snakemake.output.phases, dpi=DPI)

# Histogram of amplitudes
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.logspace(-1,4, 51)
for strictness, amplitudes_ in amplitudes.items():
    for ax, study in zip(axes.flatten(), include_studies):
        amplitude = amplitudes_[study]
        ax.hist(amplitude, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study, rotation=0, horizontalalignment="right")
        ax.set_xscale("log")
[ax.set_xlabel("Amplitude (hrs)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(snakemake.output.amplitudes, dpi=DPI)
