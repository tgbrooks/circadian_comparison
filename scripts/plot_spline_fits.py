import pandas
import numpy
import scipy.interpolate
import statsmodels.multivariate.pca
import sklearn.manifold
import pylab


# Load information
tpm = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=0).drop(columns=["Symbol"])

# Load spline fits
summary = pandas.read_csv(snakemake.input.summary, sep="\t", index_col=0)
curves = pandas.read_csv(snakemake.input.curves_fit, sep="\t", index_col=0)
curves_pstd = pandas.read_csv(snakemake.input.curves_pstd, sep="\t", index_col=0)
#re = pandas.read_csv("results/Liver/spline_fit/re.txt", sep="\t", index_col=0)
#re_structure = pandas.read_csv("results/Liver/spline_fit/re_structure.txt", sep="\t", index_col=0)

num_zeros = (tpm == 0).sum(axis=1)
summary['median_t'] =  (curves.abs()/curves_pstd).median(axis=1)

# Select genes to use
use = (summary.re_logAmp_sd < 2.5) & (summary.index.map(num_zeros) < 200) & (summary.funcDf < 10) & (summary.median_t > 2)

## Plot the gene fits for core clock genes
def alogit(x):
    return numpy.exp(x) / (numpy.exp(x) + 1)
genes = ['ENSMUSG00000055116']
names = ['Arntl']
re_by_studygene = re.reset_index().set_index(['gene', 'study'])
studies = re.study.unique()
for gene, name in zip(genes, names):
    fig, ax = pylab.subplots(figsize=(6,3))

    u = curves.loc[gene].index
    value = curves.loc[gene] + summary.loc[gene].fit_mesor

    # Plot each fit curve for each study
    for study in studies:
        logAmp, phi, mesor = re_by_genestudy.loc[(gene, study)]
        study_u = u + (alogit(phi) - 0.5)
        study_values = numpy.exp(logAmp)*value + mesor
        fig.plot((study_u%1)*24, study_values, color='gray')

    # Plot the overall fit
    fig.plot(u*24, value, color='k')

    fig.title(f"{gene} | {name}")
    fig.savefig(snakemake.output.gene_plot_dir+f"/{gene}.png")


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
yscale = (numpy.max(pca.scores['comp_1']) - numpy.min(pca.scores['comp_1'])) * 0.03
colormap = pylab.get_cmap('viridis')
colorby = 'median_t'
scale = numpy.log10
cmin, cmax = scale(summary.loc[use, colorby].min()), scale(summary.loc[use, colorby].max())
for gene, curve in curves_normalized.sample(300).iterrows():
    h = ax.plot(
             xscale * (curve.index - 0.5) + pca.scores.loc[gene, 'comp_0'],
             yscale * curve.values + pca.scores.loc[gene, 'comp_1'],
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
yscale = (numpy.max(tsne['comp_1']) - numpy.min(tsne['comp_1'])) * 0.05
colormap = pylab.get_cmap('viridis')
colorby = 'median_t'
scale = numpy.log10
cmin, cmax = scale(summary.loc[use, colorby].min()), scale(summary.loc[use, colorby].max())
for gene, curve in curves_normalized.sample(300).iterrows():
    h = ax.plot(
             xscale * (curve.index - 0.5) + tsne.loc[gene, 'comp_0'],
             yscale * curve.values + tsne.loc[gene, 'comp_1'],
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

