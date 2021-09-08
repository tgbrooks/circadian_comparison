'''
The consensus PCA is a PCA meta-analysis computed by finding the subspaces that
capture the largest variance *within* groupings, rather than between groups.
We perform this on our data, grouping by study.

Consensus PCA works by normalizing within each group before extracting 
the singular vectors, thereby finding only 'within-group' variance.
'''
import pathlib
import math

import pandas
import numpy
import scipy.sparse.linalg
import pylab

tissue = snakemake.wildcards.tissue

outdir = pathlib.Path(f"results/{tissue}/consensus_pca/")
outdir.mkdir(exist_ok=True)

data = pandas.read_csv(f"results/{tissue}/tpm_all_samples.txt", sep="\t", index_col=[0,1])
data = numpy.log(data + 0.01) # Log transform

sample_info = pandas.read_csv(f"results/{tissue}/all_samples_info.txt", sep="\t", index_col=0)
assert (data.columns == sample_info.index).all()

outlier_samples = [x.strip() for x in open(f"results/{tissue}/outlier_samples.txt").readlines()]
print(f"Dropping {len(outlier_samples)} outlier samples")

data = data.drop(columns=outlier_samples)
sample_info = sample_info.drop(index=outlier_samples)

robustness = pandas.read_csv(f"results/{tissue}/robustness_score.txt", sep="\t", index_col=0)
assert (robustness.index == data.index.get_level_values(0)).all()

#robust_genes = robustness >= 0# (sample_info.study.nunique()) * 0.5
#robust_data = data.loc[robust_genes.values]



def standardize(x):
    ''' z-score data '''
    std =  x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
    std[std.isna().any(axis=1)] = 0 # Constant rows should give 0s not NaNs
    std[x.std(axis=1) < 1e-10] = 0
    return std

def normalize_singularvalues(x):
    ''' scale top singularvalue to 1 '''
    s = numpy.linalg.svd(x.values, compute_uv=False)
    if s[0] == 0:
        return x
    return x / s[0]

# Standardize within groups
std_data = pandas.concat(
        [standardize(x) for study,x in data.groupby(sample_info.study, axis=1)],
        axis=1)

# Scale so that all studies have equal weight
normalized_data = pandas.concat(
        [normalize_singularvalues(x) for study, x in std_data.groupby(sample_info.study, axis=1)],
        axis=1)

# Perform the consensus PCA
u, s, vt = scipy.sparse.linalg.svds(normalized_data, k=4)

# Plot all the projections

nstudies = sample_info.study.nunique()
ncols = 8
nrows = math.ceil(nstudies / ncols)
remove_unplotted = slice(0,0) if nstudies == ncols * nrows else slice(-(nrows* ncols) + nstudies, None)
pca_components = {
    ("First", "Second"): (-1,-2),
    ("Second", "Third"): (-2, -3),
    ("Third", "Fourth"):(-3,- 4)
}
for (labelA, labelB), (A,B) in pca_components.items():
    fig, axes = pylab.subplots(
        figsize=(2.5*nrows, 2.5*nrows),
        ncols=ncols, nrows=nrows,
        sharex=True, sharey=True
        )
    cmap = pylab.get_cmap("twilight")
    for (study, ax) in zip(sample_info.study.unique(), axes.flatten()):
        scores = u.T @ std_data.loc[:, std_data.columns.map(sample_info.study) == study]
        times = sample_info[sample_info.study == study].time % 24
        ax.scatter(scores.values[A], scores.values[B], c=cmap(times / 24))
        ax.set_title(study)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.flatten()[remove_unplotted]:
        ax.remove()
    fig.savefig(outdir / f"consensus_pca.{labelA}.{labelB}.png", dpi=300)

# Plot the gene loadings
fig, ax = pylab.subplots(figsize=(6,6))
cmap = pylab.get_cmap("viridis")
h = ax.scatter(
    u[:,-2],
    u[:,-3],
    #c =  cmap(robustness['0'] / max(robustness)))
    c = robustness['0'])
#fig.colorbar(h)
fig.savefig(outdir / f"consensus_pca.gene_weights.Second.Third.png", dpi=300)

# Save gene loadings to disk
loadings = pandas.DataFrame({
    "1": u[:,-1],
    "2": u[:,-2],
    "3": u[:,-3],
    "4": u[:,-4],
}, index=data.index)
loadings.to_csv(outdir / "consensus_loadings.txt", sep="\t")


# Compute and plot PCAs of each study individually
# I.e. instead of a consesus PCA, we just perform PCA on 
# each study, ignoring the rest.
for (labelA, labelB), (A,B) in pca_components.items():
    fig, axes = pylab.subplots(
        figsize=(2.5*nrows, 2.5*nrows),
        ncols=ncols, nrows=nrows,
        sharex=True, sharey=True
        )
    cmap = pylab.get_cmap("twilight")
    for (study, ax) in zip(sample_info.study.unique(), axes.flatten()):
        study_data = data.loc[:,sample_info.study == study]
        std_study_data = standardize(study_data)
        u2, s2, vt2 = scipy.sparse.linalg.svds(std_study_data, k=4)
        scores = u2.T @ std_study_data
        times = sample_info[sample_info.study == study].time % 24
        ax.scatter(scores.values[A], scores.values[B], c=cmap(times / 24))
        ax.set_title(study)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.flatten()[remove_unplotted]:
        ax.remove()
    fig.savefig(outdir / f"individual_pca.{labelA}.{labelB}.png", dpi=300)
