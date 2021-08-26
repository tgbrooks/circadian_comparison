import pathlib
import pandas
import numpy
import scipy.interpolate
import statsmodels.multivariate.pca
import sklearn.manifold
import pylab

DPI = 300

# Load information
#tpm = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=0).drop(columns=["Symbol"])

# Load spline fits
summary = pandas.read_csv(snakemake.input.summary, sep="\t", index_col=0)
curves = pandas.read_csv(snakemake.input.curves_fit, sep="\t", index_col=0)
curves_pstd = pandas.read_csv(snakemake.input.curves_pstd, sep="\t", index_col=0)
re = pandas.read_csv("results/Liver/spline_fit/re.txt", sep="\t", index_col=0) # random effects
re_structure = pandas.read_csv("results/Liver/spline_fit/re_structure.txt", sep="\t", index_col=0)

#num_zeros = (tpm == 0).sum(axis=1)
summary['median_t'] =  (curves.abs()/curves_pstd).median(axis=1)
summary['rel_amp'] = summary['fit_amplitude'] / summary['sigma']

# Select genes to use
use = (summary.re_logAmp_sd < 3) & (summary.funcDf < 15) & (summary.median_t > 2) & (summary.iteration_times < 31)
print(f"Using {use.sum()} out of {len(use)} genes.")

## Plot the gene fits for core clock genes
def alogit(x):
    return numpy.exp(x) / (numpy.exp(x) + 1)
genes = ["ENSMUSG00000055116", "ENSMUSG00000020038", "ENSMUSG00000068742", "ENSMUSG00000020893", "ENSMUSG00000055866", "ENSMUSG00000028957", "ENSMUSG00000020889", "ENSMUSG00000021775", "ENSMUSG00000059824", "ENSMUSG00000029238"]
names= ["Arntl", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Nr1d1", "Nr1d2", "Dbp", "Clock"]
re_by_studygene = re.reset_index().set_index(['gene', 'study'])
studies = re.study.unique()
pathlib.Path(snakemake.output.gene_plot_dir).mkdir(exist_ok=True)
for gene, name in zip(genes, names):
    if gene not in summary.index:
        print(f"Gene {gene} | {name} has no spline fit.")
        continue
    fig, ax = pylab.subplots(figsize=(6,3))

    u = curves.loc[gene].index.astype(float)
    value = curves.loc[gene] + summary.loc[gene].fit_mesor

    # Plot each fit curve for each study
    for study in studies:
        logAmp, phi, mesor = re_by_studygene.loc[(gene, study)]
        study_u = ((u + (alogit(phi) - 0.5)) % 1)
        study_values = numpy.exp(logAmp)*value + mesor
        order = numpy.argsort(study_u)
        ax.plot( (study_u*24)[order], numpy.exp(study_values[order]), color='gray', linewidth=0.5)

    # Plot the overall fit
    ax.plot(u*24, numpy.exp(value), color='k')

    ax.set_title(f"{gene} | {name}")
    fig.savefig(snakemake.output.gene_plot_dir+f"/{gene}.png", dpi=DPI)

# Plot the phase distributions
fig, axes = pylab.subplots(figsize=(5,5,), nrows=2, sharex=True, sharey=True)
bins = numpy.linspace(0,24,25) # hourly bins
axes[0].hist(summary[use].fit_peak_time)
axes[0].set_ylabel("Peak Time")
axes[1].hist(summary[use].fit_trough_time)
axes[1].set_ylabel("Trough Time")
axes[1].set_xlabel("Time (hours)")
axes[1].set_xticks(numpy.linspace(0,24,5))
fig.tight_layout()
fig.savefig(snakemake.output.phase_distribution, dpi=DPI)

# Plot correlation of trough and peak time
fig, ax = pylab.subplots()
ax.scatter(
        summary[use].fit_peak_time + 0.1*numpy.random.normal(size=len(summary[use])),
        summary[use].fit_trough_time + 0.1*numpy.random.normal(size=len(summary[use])),
        s=1)
ax.plot([0,12],[12,24], c='k', linewidth=1, zorder=-1)
ax.plot([12,24],[0,12], c='k', linewidth=1, zorder=-1)
ax.set_xlim(0,24)
ax.set_ylim(0,24)
ax.set_xlabel("Peak Time (hrs)")
ax.set_ylabel("Trough Time (hrs)")
ax.set_xticks(numpy.linspace(0,24,5))
ax.set_yticks(numpy.linspace(0,24,5))
fig.savefig(snakemake.output.phase_correlation, dpi=DPI)

# Align the curves by their peak time
# So all peaks will occur at time u=0.5 (out of the range [0,1])
curves_long = curves[use].stack().rename_axis(index=["gene", "u"]).reset_index()
curves_long.u = curves_long.u.astype(float)
curves_long_shifted = curves_long.copy() # Shift u according to the peak time
curves_long_shifted.u = (curves_long_shifted.u -  curves_long_shifted.gene.map(summary.fit_peak_time) /24 + 0.5) % 1.0
aligning = []
for gene, data in curves_long_shifted.groupby("gene"):
    aligned =  pandas.DataFrame({
        "gene": gene,
        "u": numpy.linspace(0,1, 24*4+1),
    })
    # triple the data cyclically so that we can always interpolate
    data1 = data.copy()
    data1.u -= 1
    data2 = data.copy()
    data2.u += 1
    data = pandas.concat([
        data1,
        data,
        data2
    ])
    aligned['value'] = scipy.interpolate.interp1d(data.u, data[0], fill_value="extrapolate")(aligned.u)
    aligning.append(aligned)
curves_aligned = pandas.concat(aligning)
curves_aligned = curves_aligned.pivot(index="gene", columns="u", values="value")

# Scale to be [-1,1]
curves_normalized = curves_aligned.div(curves_aligned.abs().max(axis=1), axis=0)

# Temp fix to drop crazy outliers
outliers = ['ENSMUSG00000064356']
curves_normalized = curves_normalized[~curves_normalized.index.isin(outliers)]

# Create a PCA plot of the curves
pca = statsmodels.multivariate.pca.PCA(curves_normalized, ncomp=2)


fig, ax = pylab.subplots(figsize=(8,8))
#ax.scatter(pca.scores['comp_0'], pca.scores['comp_1'])
xscale = (numpy.max(pca.scores['comp_0']) - numpy.min(pca.scores['comp_0'])) * 0.08
yscale = (numpy.max(pca.scores['comp_1']) - numpy.min(pca.scores['comp_1'])) * 0.015
colormap = pylab.get_cmap('viridis')
colorby = 'rel_amp'
scale = numpy.log10
cmin, cmax = scale(summary.loc[use, colorby].min()), scale(summary.loc[use, colorby].max())
for gene, curve in curves_normalized.sample(300).iterrows():
    h = ax.plot(
             xscale * (curve.index - 0.5) + pca.scores.loc[gene, 'comp_0'],
             yscale * numpy.exp(curve.values) + pca.scores.loc[gene, 'comp_1'],
             c=colormap((scale(summary.loc[gene, colorby]) - cmin) / (cmax - cmin)),
             linewidth=0.5)
ax.set_axis_off()
ax.set_title("PCA of Curve Fits")
#fig.colorbar(h)
fig.savefig(snakemake.output.pca, dpi=300)

# Run a tsne dimension reductin on the splines
tsne = sklearn.manifold.TSNE().fit_transform(curves_normalized)
tsne = pandas.DataFrame(tsne, index=curves_normalized.index, columns=['comp_0', 'comp_1'])

fig, ax = pylab.subplots(figsize=(8,8))
#ax.scatter(tsne['comp_0'], tsne['comp_1'])
# Downsample the genes and plot their curves
xscale = (numpy.max(tsne['comp_0']) - numpy.min(tsne['comp_0'])) * 0.05
yscale = (numpy.max(tsne['comp_1']) - numpy.min(tsne['comp_1'])) * 0.025
colormap = pylab.get_cmap('viridis')
colorby = 'rel_amp'
scale = numpy.log10
cmin, cmax = scale(summary.loc[use, colorby].min()), scale(summary.loc[use, colorby].max())
for gene, curve in curves_normalized.sample(300).iterrows():
    h = ax.plot(
             xscale * (curve.index - 0.5) + tsne.loc[gene, 'comp_0'],
             yscale * numpy.exp(curve.values) + tsne.loc[gene, 'comp_1'],
             c=colormap((scale(summary.loc[gene, colorby]) - cmin) / (cmax - cmin)),
             linewidth=0.5)
ax.set_axis_off()
ax.set_title("t-SNE of Curve Fits")
fig.tight_layout()
#fig.colorbar(h)
fig.savefig(snakemake.output.tsne, dpi=300)


## Just plot the curves ontop of eachother
# maybe not so useful
#fig, ax = pylab.subplots(figsize=(3,20))
#for idx, (gene, curve) in enumerate(curves_normalized.sample(500).iterrows()):
#    ax.plot( curve.index, idx + curve.values, c='k', linewidth=0.5)
#fig.savefig("results/Liver/spline_fit/curves.png", dpi=300)

