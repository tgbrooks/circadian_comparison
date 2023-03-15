import pathlib
import math
import pandas
import numpy
import seaborn as sns
import scipy.interpolate
import scipy.special
import statsmodels.multivariate.pca
import sklearn.manifold
import pylab
import matplotlib

from studies import targets
import styles

def legend_from_colormap(fig, colormap, names=None, **kwargs):
    if names is None:
        names = {cat:cat for cat in colormap.keys()}
    legend_elts = [matplotlib.lines.Line2D(
                            [0],[0],
                            marker="o", markerfacecolor=c, markersize=10,
                            label=names[cat] if not pandas.isna(cat) else "NA",
                            c=c, lw=0)
                        for cat, c in colormap.items()]
    fig.legend(handles=legend_elts, **kwargs)

DPI = 300

# Load original information
tpm = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=0).drop(columns=["Symbol"])
sample_info = pandas.read_csv(snakemake.input.sample_info, sep="\t", index_col=0)
outlier_samples = [x.strip() for x in open(snakemake.input.outliers).readlines()]
sample_info = sample_info[~sample_info.index.isin(outlier_samples)]

# Load spline fits
summary = pandas.read_csv(snakemake.input.summary, sep="\t", index_col=0)
curves = pandas.read_csv(snakemake.input.curves_fit, sep="\t", index_col=0)
curves_pstd = pandas.read_csv(snakemake.input.curves_pstd, sep="\t", index_col=0)
re = pandas.read_csv("results/Liver/spline_fit/re.txt", sep="\t", index_col=0) # random effects
re_structure = pandas.read_csv("results/Liver/spline_fit/re_structure.txt", sep="\t", index_col=0)


statsdir = pathlib.Path(snakemake.input.statsdir)
num_peaks = pandas.read_csv(statsdir / "num_peaks.txt", sep="\t", index_col=0)['0']
asymmetric = pandas.read_csv(statsdir / "asymmetric.txt", sep="\t", index_col=0)['0']
goodness_of_fit = pandas.read_csv(statsdir / "goodness_of_fit.txt", sep="\t")

use = summary.is_rhythmic
dark2 = pylab.get_cmap("Dark2")
color_by_category = {
        "symmetric": dark2(0),
        "asymmetric": dark2(1),
        "multimodal": dark2(2),
}

## Plot the gene fits for select genes
alogit = scipy.special.expit
genes = ["ENSMUSG00000055116", "ENSMUSG00000020038", "ENSMUSG00000068742", "ENSMUSG00000020893", "ENSMUSG00000055866", "ENSMUSG00000028957", "ENSMUSG00000020889", "ENSMUSG00000021775", "ENSMUSG00000059824", "ENSMUSG00000029238", "ENSMUSG00000057342", "ENSMUSG00000016619"]
names= ["Arntl", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Nr1d1", "Nr1d2", "Dbp", "Clock", "Sphk2", "Nup50"]
re_by_studygene = re.reset_index().set_index(['gene', 'study'])
re_studies = re.study.unique()
studies = sorted(sample_info.study.unique(), key = lambda x: targets[x]['short_name'])
pathlib.Path(snakemake.output.gene_plot_dir).mkdir(exist_ok=True)
for gene, name in zip(genes, names):
    if gene not in summary.index:
        print(f"Gene {gene} | {name} has no spline fit.")
        continue

    u = curves.loc[gene].index.astype(float)
    value = curves.loc[gene]

    # Plot the original data aligned by the fit phase/amplitude
    # and with the overall fit ontop
    ncols = 7
    nrows = math.ceil(len(studies) / ncols)
    fig, axes = pylab.subplots(figsize=(12,12), nrows=nrows, ncols=ncols, sharex=True, sharey=True)
    for study, ax in zip(studies, axes.flatten()):
        study_samples = sample_info.index[sample_info.study == study]
        study_tpm = tpm.loc[gene, study_samples]
        times = sample_info.loc[study_samples].time
        # Plot raw data
        ax.scatter(times%24, numpy.log(study_tpm+0.01), marker='+', color='r', label="Raw data")
        if study in re_studies:
            # plot the the study-specific spline fit
            logAmp, phi, mesor = re_by_studygene.loc[(gene, study)]
            study_u = ((u + (alogit(phi) - 0.5)) % 1)
            study_values = numpy.exp(logAmp)*value + mesor + summary.loc[gene].fit_mesor
            order = numpy.argsort(study_u)
            ax.plot((study_u * 24)[order], study_values[order], color='k', label='Study fit')
        ax.set_title(targets[study]['short_name'])
        ax.set_xticks([0,6,12,18,24])
    for ax in axes[:,0]:
        ax.set_ylabel("log TPM")
    for col in range(ncols):
        # Label the x axis, taking into account how there may be
        # ragged columns and so sometimes we need to label the second
        # to bottom row
        overhang = len(studies) % ncols # num cols on bottom row
        in_last_row = (overhang > col) or (overhang == 0)
        if in_last_row:
            row_idx = nrows - 1
        else:
            row_idx = max(0, nrows - 2)
        axes[row_idx,col].xaxis.set_tick_params(labelbottom=True)
        axes[row_idx,col].set_xlabel("Time")
    for ax in axes.flatten()[len(studies):]:
        ax.remove()
    fig.suptitle(f"{gene} | {name}")
    fig.tight_layout()
    fig.savefig(snakemake.output.gene_plot_dir+f"/{gene}.by_study.png", dpi=DPI)
    fig.savefig(snakemake.output.gene_plot_dir+f"/{gene}.by_study.svg")
    pylab.close(fig)

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
pylab.close(fig)

# Plot correlation of trough and peak time
joint = sns.jointplot(x="fit_peak_time", y="fit_trough_time", data=summary.loc[use], kind="hex")
joint.ax_joint.plot([0,12],[12,24], c='k', linewidth=1)
joint.ax_joint.plot([12,24],[0,12], c='k', linewidth=1)
joint.ax_joint.set_xlabel("Acrophase (hrs)")
joint.ax_joint.set_ylabel("Bathyphase (hrs)")
joint.ax_joint.set_xticks(numpy.linspace(0,24,5))
joint.ax_joint.set_yticks(numpy.linspace(0,24,5))
joint.figure.savefig(snakemake.output.phase_correlation, dpi=DPI)
pylab.close(joint.figure)

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
             yscale * curve.values + pca.scores.loc[gene, 'comp_1'],
             #c=colormap((scale(summary.loc[gene, colorby]) - cmin) / (cmax - cmin)),
             c = color_by_category[summary.loc[gene, 'category']],
             linewidth=1.5)
ax.set_axis_off()
ax.set_title("PCA of Curve Fits")
legend_from_colormap(fig, color_by_category)
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
             yscale * (curve.values) + tsne.loc[gene, 'comp_1'],
             #c=colormap((scale(summary.loc[gene, colorby]) - cmin) / (cmax - cmin)),
             c = color_by_category[summary.loc[gene, 'category']],
             linewidth=1.5)
ax.set_axis_off()
ax.set_title("t-SNE of Curve Fits")
fig.tight_layout()
legend_from_colormap(fig, color_by_category)
#fig.colorbar(h)
fig.savefig(snakemake.output.tsne, dpi=300)
