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

# Load robustness scores and TPM
robustness_score = pandas.read_csv(snakemake.input.robustness, sep="\t", index_col=0)["0"]
data = {}
for study, tpmfile in zip(studies, snakemake.input.tpm):
    tpm = pandas.read_csv(tpmfile, sep="\t", index_col=0)
    data[study] = tpm

data_df = pandas.concat(data.values(), axis=1)

med_tpm = data_df.median(axis=1)

# Load JTK vlaues
amp = {}
per = {}
phase = {}
p = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    jtk.loc[(jtk['ADJ.P'] >= 0.05) | (jtk['dropped']), 'LAG'] = float("nan")
    jtk.loc[(jtk['ADJ.P'] >= 0.05) | (jtk['dropped']), 'PER'] = float("nan")
    amp[study] = jtk['AMP']
    per[study] = jtk["PER"]
    p[study] = jtk["ADJ.P"]
    phase[study] =jtk["LAG"]

amp_df = pandas.concat(amp.values(), axis=1)
per_df = pandas.concat(per.values(), axis=1)
phase_df = pandas.concat(phase.values(), axis=1)
p_df = pandas.concat(p.values(), axis=1)

# Robustness groups
robustness_cat = pandas.cut(
    robustness_score,
    bins = numpy.arange(0,robustness_score.max()+5,5),
    include_lowest=True,
)

# Robustness score histogram
fig, ax = pylab.subplots(figsize=(4,4))
ax.hist(robustness_score)
ax.set_xlabel("Robustness Score")
ax.set_ylabel("Num. Genes")
ax.set_yscale('log')
fig.tight_layout()
fig.savefig(snakemake.output.histogram, dpi=DPI)


# Robustness score versus expression level (TPM)
fig, ax = pylab.subplots(figsize=(4,4))
ax.scatter(
    med_tpm,
    robustness_score+numpy.random.normal(size=len(robustness_score))*0.2,
    s=1,
    alpha=0.3)
ax.set_xscale("log")
ax.set_xlabel("Median TPM")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.expression_level, dpi=DPI)

#amp, period, phase from jtk

# robustness verus amplitude
med_amp = amp_df.median(axis=1)
fig, ax = pylab.subplots(figsize=(4,4))
ax.scatter(
    med_amp,
    robustness_score+numpy.random.normal(size=len(robustness_score))*0.2,
    s=1,
    alpha=0.3)
ax.set_xscale("log")
ax.set_xlabel("Median Amplitude")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.amplitude, dpi=DPI)

# Robustness versus period
med_per = per_df.median(axis=1)
fig, ax = pylab.subplots(figsize=(4,4))
score_by_period = dict(list(robustness_score.groupby(med_per)))
period_by_score = dict(list(med_per.groupby(robustness_cat)))
bins = numpy.linspace(18.5,24.5,7)
ax.scatter(
    med_per + numpy.random.uniform(-0.4, 0.4, size=len(med_per)),
    robustness_score,
    alpha=0.3,
    s=1,
)
#for score_cat, periods in period_by_score.items():
#    ax.hist(
#        periods,
#        bins = bins,
#        #density=True,
#        weights=[0.9*(score_cat.right - score_cat.left) / periods.count()
#                    for i in range(len(periods))],
#        bottom = score_cat.left,
#        color='k',
#    )
#parts = ax.violinplot(
#    #list(score_by_period.values()),
#    #positions=list(score_by_period.keys()),
#    list(period_by_score.values()),
#    positions=[interval.mid for interval in period_by_score.keys()],
#    showextrema=False,
#    widths=2.0,
#    vert=False,
#)
#for body in parts['bodies']:
#    body.set_alpha(1)
ax.set_xlim(19,25)
ax.set_xlabel("Median Period")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.period, dpi=DPI)

# By phase
def round(x, to=1):
    return numpy.round(x/to)*to
mean_phase = scipy.stats.circmean(phase_df, low=0, high=24, axis=1, nan_policy="omit")
mean_phase = pandas.Series(mean_phase, index=robustness_score.index)
fig, ax = pylab.subplots(figsize=(5,4))
#score_by_phase = dict(list(robustness_score.groupby(round(mean_phase, to=2) % 24)))
#phase_by_score = dict(list(mean_phase.groupby(robustness_cat)))
#parts = ax.violinplot(
#    #list(score_by_phase.values()),
#    #positions=list(score_by_phase.keys()),
#    list(phase_by_score.values()),
#    positions=[interval.mid for interval in phase_by_score.keys()],
#    showextrema=False,
#    widths=2.0,
#    vert=False,
#)
#for body in parts['bodies']:
#    body.set_alpha(0.5)
ax.scatter(mean_phase, robustness_score, s=1, alpha=0.3, zorder=-5)
ax.set_xticks([0,6,12,18,24])
ax.set_xlabel("Mean Phase")
ax.set_ylabel("Robustness Score")
fig.savefig(snakemake.output.phase, dpi=DPI)

