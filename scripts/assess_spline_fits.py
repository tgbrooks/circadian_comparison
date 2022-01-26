import pathlib
import math
import pandas
import numpy
import scipy.interpolate
import statsmodels.multivariate.pca
import sklearn.manifold
import pylab

import styles

DPI = 300

# Cutoff point at which we consider the t-stat of the
# fit curves to be elevated above the mean
T_STAT_ELEVATED_CUTOFF = 2

out_dir = pathlib.Path(snakemake.output[0])
out_dir.mkdir(exist_ok=True)

# Load original information
tpm = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=0).drop(columns=["Symbol"])
sample_info = pandas.read_csv(snakemake.input.sample_info, sep="\t", index_col=0)
outlier_samples = [x.strip() for x in open(snakemake.input.outliers).readlines()]
sample_info = sample_info[~sample_info.index.isin(outlier_samples)]

# Load spline fits
summary = pandas.read_csv(snakemake.input.summary, sep="\t", index_col=0)
curves = pandas.read_csv(snakemake.input.curves_fit, sep="\t", index_col=0)
curves_pstd = pandas.read_csv(snakemake.input.curves_pstd, sep="\t", index_col=0)
re = pandas.read_csv(snakemake.input.re, sep="\t", index_col=0) # random effects
re_structure = pandas.read_csv(snakemake.input.re_structure, sep="\t", index_col=0)

#num_zeros = (tpm == 0).sum(axis=1)
summary['median_t'] =  (curves.abs()/curves_pstd).median(axis=1)
summary['rel_amp'] = summary['fit_amplitude'] / summary['sigma']

# Select genes to use
use = (summary.re_logAmp_sd < 3) & (summary.funcDf < 15) & (summary.median_t > 2) & (summary.iteration_times < 31)
print(f"Using {use.sum()} out of {len(use)} genes.")


# Setup
re_by_studygene = re.reset_index().set_index(['gene', 'study'])
re_studies = re.study.unique()
studies = [study for study in styles.studies if study in re_studies] # Get order of studies consistent

def alogit(x):
    return numpy.exp(x) / (numpy.exp(x) + 1)

#### Try to classify the curves into different categories
# First: multi-modal or mono-modal
# We will define multi-modal by looking at the curves / curves_pstd 't-scores'
# and looking for two peaks above mean, or two peaks below mean
# Specifically we need to find two points that are both significantly elevated above
# the mean (0) and which have a value < 0 between them (not necessarily significantly so)
tstats = curves / curves_pstd
def count_peaks(tstats):
    def count_peaks_row(row):
        first_peak = None
        num_switches = 0 # Number of times we switch from peak to trough
        next_target = "peak" # Whether we need a peak or trough next
        i = 0
        # Iterate through a full cycle, starting at the first identified peak (if any)
        while True:
            if row.iloc[i % len(row)] > T_STAT_ELEVATED_CUTOFF:
                if first_peak is None:
                    first_peak = i
                if next_target == 'peak':
                    next_target = 'trough'
            if row.iloc[i % len(row)] < 0:
                if next_target == 'trough':
                    next_target = 'peak'
                    num_switches += 1
            i += 1
            if first_peak is None:
                if i >= len(row):
                    break
            else:
                if i >= len(row) + first_peak:
                    break
        return num_switches

    num_peaks = tstats.apply(count_peaks_row, axis=1)
    num_troughs = (-tstats).apply(count_peaks_row, axis=1)

    return numpy.maximum(num_peaks, num_troughs)

num_peaks = count_peaks(tstats)
num_peaks.to_csv(out_dir / "num_peaks.txt", sep="\t")

def calc_goodness_of_fit(curves, tpm):
    res = []
    for gene in curves.index:
        u = curves.loc[gene].index.astype(float)
        value = curves.loc[gene]

        for study in studies:
            logAmp, phi, mesor = re_by_studygene.loc[(gene, study)]
            study_samples = sample_info.index[sample_info.study == study]
            study_tpm = tpm[study_samples].loc[gene]
            study_u = (sample_info.loc[study_samples].time / 24 - (alogit(phi) - 0.5))
            study_values = numpy.log(study_tpm+0.01)
            fit_values = numpy.exp(logAmp)*numpy.interp(study_u, u, value, period=1) + mesor + summary.loc[gene].fit_mesor
            resid = (study_values - fit_values)
            res.append({
                "gene": gene,
                "study": study,
                "SS_resid": (resid**2).sum(),
                "SS_total": ((study_values - study_values.mean())**2).sum(),
                "SS_reg": ((fit_values - fit_values.mean())**2).sum(),
            })
    res = pandas.DataFrame(res)
    res['Rsquared'] = res.SS_reg / res.SS_total
    res['Rsquared2'] = (res.SS_total - res.SS_resid) / (res.SS_total)
    return res
goodness_of_fit = calc_goodness_of_fit(curves, tpm)
goodness_of_fit.to_csv(out_dir / "goodness_of_fit.txt", sep="\t")
