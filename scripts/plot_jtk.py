import json
import math
import pathlib
import pandas
import numpy
import scipy.stats
import matplotlib
matplotlib.use("Agg")
import pylab

import util
from styles import format_study_name, color_by_sex, color_by_light
from studies import targets as study_info

# Order studies alphabetically
studies = sorted(snakemake.params.studies, key = lambda x: study_info[x]['short_name'])
N_studies = len(studies)

robustness_score = pandas.read_csv(snakemake.input.robustness_score, sep="\t", index_col=0)['0']
highly_robust_genes = robustness_score.index[robustness_score >= 35]

jtk_period = snakemake.wildcards.period
if jtk_period == '':
    jtk_period = 24
else:
    jtk_period = int(jtk_period)

DPI = 300
N_COLUMNS = 3
# Q-value cutoff to call significantly rhythmic
Q_CUTOFFS = {
    "loose": 0.05,
    "strict": 0.01,
}
color_by_strictness = {
    "loose": "lightblue",
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
robust_genes_phases = {}
amplitudes = {strictness:{} for strictness in Q_CUTOFFS.keys()}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    breakdowns[study] = num_below(jtk['qvalue'], breakpoints)

    for strictness, cutoff in Q_CUTOFFS.items():
        significant = jtk[jtk['qvalue'] < cutoff]
        periods[strictness][study] = significant['PER']
        phases[strictness][study] = significant['LAG']
        amplitudes[strictness][study] = significant['AMP']

        robust_genes_phases[study] = jtk.loc[jtk.index.isin(highly_robust_genes), 'LAG']


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
if len(include_studies) == 0:
    # none passed threshold, might as well plot all
    include_studies = studies
N_studies = len(include_studies)
N_ROWS = math.ceil(N_studies / N_COLUMNS)
remove_unplotted = slice(0,0) if N_studies == N_ROWS * N_COLUMNS else slice(-(N_ROWS * N_COLUMNS) + N_studies, None)

# Histogram of periods
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.linspace(19.5,28.5,10)
for strictness, periods_ in periods.items():
    for ax, study in zip(axes.flatten(), include_studies):
        period = periods_[study]
        ax.hist(period, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study_info[study]['short_name'], rotation=0, horizontalalignment="right")
[ax.set_xlabel("Period (hrs)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(snakemake.output.periods, dpi=DPI)

# Histogram of phases
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.linspace(0,jtk_period,jtk_period+1)
for strictness, phases_ in phases.items():
    for ax, study in zip(axes.flatten(), include_studies):
        phase = phases_[study]
        ax.hist(phase, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study_info[study]['short_name'], rotation=0, horizontalalignment="right")
        ax.set_xlim(0,jtk_period)
        ax.set_xticks(numpy.arange(0,jtk_period+1, jtk_period//4))
[ax.set_xlabel("Phase (hrs)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(snakemake.output.phases, dpi=DPI)

# Heatmap of phases
phase_types = {
    'loose': phases['loose'],
    'robust': robust_genes_phases,
}
for phase_type, phases_of_type in phase_types.items():
    plot_studies = sorted(studies, key=lambda x: study_info[x]['short_name'])
    fig, main_ax = pylab.subplots(
            figsize=(6,9),
            constrained_layout=True
    )
    density_by_study = numpy.zeros(shape=(len(plot_studies), len(bins)-1))
    num_sig_by_study = numpy.zeros((len(plot_studies),1))
    for i, study in enumerate(plot_studies):
        phase = phases_of_type[study]
        num_sig_by_study[i,0] = len(phase)
        if len(phase) < 10:
            continue
        kde = scipy.stats.gaussian_kde(phase.values, bw_method=0.3)
        vals = kde(bins[:-1]+0.5)
        vals /= max(vals) # Maximum is 1
        density_by_study[i,:] = vals
    cmap = matplotlib.cm.get_cmap().copy()
    # main values
    h = main_ax.imshow(
        density_by_study,
        cmap = cmap,
    )
    main_ax.set_xlim(-0.5, len(bins)-1-0.5)
    time_labels = numpy.arange(0, jtk_period+1, jtk_period//4)
    main_ax.set_xticks(time_labels-0.5)
    main_ax.set_xticklabels([str(x) for x in time_labels])
    main_ax.set_xlabel("Phase (hrs)")
    main_ax.set_yticks(numpy.arange(len(plot_studies)))
    main_ax.set_yticklabels([study_info[study]['short_name'] for study in plot_studies])

    #Colorbars + legend
    fig.colorbar(h, ax=main_ax, label="Phase Density", fraction=0.025)

    if phase_types == 'loose':
        # Number of rhythmic genes labels
        dupe_ax = ax.secondary_yaxis("right")
        dupe_ax.set_yticks(numpy.arange(len(plot_studies)))
        dupe_ax.set_yticklabels(num_sig_by_study)
        dupe_ax.set_ylabel("Num. Rhythmic")

    fig.savefig(f"results/{snakemake.wildcards.tissue}/jtk{snakemake.wildcards.period}/phases.heatmap.{phase_type}.png", dpi=DPI)
    fig.savefig(f"results/{snakemake.wildcards.tissue}/jtk{snakemake.wildcards.period}/phases.heatmap.{phase_type}.svg")

# Histogram of amplitudes
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.logspace(-1,4, 51)
for strictness, amplitudes_ in amplitudes.items():
    for ax, study in zip(axes.flatten(), include_studies):
        amplitude = amplitudes_[study]
        ax.hist(amplitude, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study_info[study]['short_name'], rotation=0, horizontalalignment="right")
        ax.set_xscale("log")
[ax.set_xlabel("Amplitude (TPM)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(snakemake.output.amplitudes, dpi=DPI)
