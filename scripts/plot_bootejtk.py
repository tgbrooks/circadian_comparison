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

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

robustness_score = pandas.read_csv(snakemake.input.robustness_score, sep="\t", index_col=0)['0']
highly_robust_genes = robustness_score.index[robustness_score >= 35]

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
jtk_period = 24

# Load bootejtk data
bootejtk = pandas.read_csv(snakemake.input.bootejtk, sep="\t")

# Order studies alphabetically
studies = sorted(bootejtk.study.unique(), key = lambda x: study_info[x]['short_name'])
N_studies = len(studies)

N_ROWS = math.ceil(N_studies / N_COLUMNS)
remove_unplotted = slice(0,0) if N_studies == N_ROWS * N_COLUMNS else slice(-(N_ROWS * N_COLUMNS) + N_studies, None)

# Histogram of phases
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.linspace(0,jtk_period,jtk_period+1)
for strictness, q_cutoff in Q_CUTOFFS.items():
    for (study, study_data), ax in zip(bootejtk.groupby("study"), axes.flatten()):
        phase = study_data[study_data.GammaBH < q_cutoff].PhaseMean
        ax.hist(phase, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study_info[study]['short_name'], rotation=0, horizontalalignment="right")
        ax.set_xlim(0,jtk_period)
        ax.set_xticks(numpy.arange(0,jtk_period+1, jtk_period//4))
[ax.set_xlabel("Phase (hrs)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(outdir / "phases.png", dpi=DPI)

# Heatmap of phases
phase_types = {
    'loose': bootejtk[bootejtk.GammaBH < Q_CUTOFFS['loose']][['study', 'ID', 'PhaseMean']],
    'robust': bootejtk[bootejtk.ID.isin(highly_robust_genes)][['study', 'ID', 'PhaseMean']],
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
        phase = phases_of_type[phases_of_type.study == study].PhaseMean
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

    fig.savefig(outdir / f"phases.heatmap.{phase_type}.png", dpi=DPI)
    fig.savefig(outdir / f"phases.heatmap.{phase_type}.svg")

# Histogram of amplitudes
fig, axes = pylab.subplots(figsize=(1+5*N_COLUMNS,0.7+0.7*N_ROWS), nrows=N_ROWS, ncols=N_COLUMNS, sharex=True, squeeze=False)
bins = numpy.logspace(-1,4, 51)
for strictness, q_cutoff in Q_CUTOFFS.items():
    for (study, study_data), ax in zip(bootejtk.groupby("study"), axes.flatten()):
        amplitude = study_data[study_data.GammaBH < q_cutoff].Max_Amp
        ax.hist(amplitude, bins=bins, color=color_by_strictness[strictness])
        ax.set_ylabel(study_info[study]['short_name'], rotation=0, horizontalalignment="right")
        ax.set_xscale("log")
[ax.set_xlabel("Amplitude (read counts)") for ax in axes[-1,:]]
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
util.legend_from_colormap(fig, color_by_strictness, names={s:f"Q < {c:0.2f}" for s,c in Q_CUTOFFS.items()})
fig.tight_layout()
fig.savefig(outdir / "amplitudes.png", dpi=DPI)
