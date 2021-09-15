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
#tissue = "Liver"

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


# Utility functions
# Plot all the projections
def format_study(study, max_length=15):
    parts = study.split("_")
    line = parts[0]
    lines = []
    for part in parts[1:]:
        if len(line) + len(part) + 1 > max_length:
            lines.append(line)
            line = part
        else:
            line = line + "_" + part
    lines.append(line)
    return '\n'.join(lines)


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

## Perform a "JIVE" analysis (see PMID: 23745156)
def JIVE(data, groupings, RJ=2, RA=1, MAX_ITER=3):
    ''' Given normalized data, compute JIVE by alternating between
    identifying joint and individual #CAs

    data = normalized data
    groups = array of group assignments for the columns
    RJ = rank of the joint matrices J_i
    RA = rank of the individual matrices A_i'''

    print("Starting JIVE")
    # Initial loading is to have no joint factor
    joint_loading = numpy.zeros(shape=(RJ, data.shape[0]))

    for i in range(MAX_ITER):
        # Project data to the joint factors
        joint_matrix = joint_loading.T @ (joint_loading @ data)

        joint_residual = data.values - joint_matrix

        # Find the individual groups' PCA and remove those
        individual_residuals = []
        for group, group_resid in joint_residual.groupby(groupings, axis=1):
            print('.', end='', flush=True)
            # Compute PCA after removing joint components
            u, s, vt = numpy.linalg.svd(group_resid, full_matrices=False)
            individual_loadings = u[:,:RA].T #Top components

            # Compute residual without the joint
            group_data = data.loc[:, groupings == group]
            individual_matrix = individual_loadings.T @ (individual_loadings @ group_data)
            residual = group_data.values - individual_matrix
            individual_residuals.append(residual)
        individual_residual = numpy.concatenate(individual_residuals, axis=1)

        # Compute the joint PCA from the residual from individual PCAs
        u, s, vt = numpy.linalg.svd(individual_residual, full_matrices=False)
        joint_loading = u[:,:RJ].T
        total_resid = individual_residual - joint_loading.T @ (joint_loading @ individual_residual)
        print(f"Iteration {i} complete")
        print(f"total residual: {numpy.linalg.norm(total_resid):0.3f} of {numpy.linalg.norm(data):0.3f}")

    joint_singular_vals = numpy.linalg.norm(joint_loading @ data, axis=1)
    joint_percent_explained = joint_singular_vals**2 / numpy.linalg.norm(data)**2

    return pandas.DataFrame(joint_loading, columns=data.index), joint_singular_vals, joint_percent_explained

# Compute JIVE
# We use just 2 joint components and 1 individual component
joint_loading, joint_sing_vals, joint_pct_explained = JIVE(normalized_data, sample_info.study, RJ=2, RA=1)

joint_loading.to_csv(outdir/ "JIVE_joint_loadings.txt", sep="\t")
pandas.Series(joint_pct_explained).to_csv(outdir / "JIVE_joint_pct_explained.txt", sep="\t")

# Plot JIVE
nstudies = sample_info.study.nunique()
ncols = 8
nrows = math.ceil(nstudies / ncols)
remove_unplotted = slice(0,0) if nstudies == ncols * nrows else slice(-(nrows* ncols) + nstudies, None)
fig, axes = pylab.subplots(
    figsize=(2.5*nrows, 2.5*nrows),
    ncols=ncols, nrows=nrows,
    sharex=True, sharey=True
    )
cmap = pylab.get_cmap("twilight")
for (study, ax) in zip(sample_info.study.unique(), axes.flatten()):
    scores = joint_loading @ std_data.loc[:, std_data.columns.map(sample_info.study) == study]
    times = sample_info[sample_info.study == study].time % 24
    ax.scatter(scores.values[0], scores.values[1], c=cmap(times / 24))
    ax.set_title(format_study(study))
    ax.set_xticks([])
    ax.set_yticks([])
for ax in axes.flatten()[remove_unplotted]:
    ax.remove()
fig.supxlabel(f"Joint PC 1 ({joint_pct_explained[0]:0.1%})")
fig.supylabel(f"Joint PC 2 ({joint_pct_explained[1]:0.1%})")
fig.tight_layout()
fig.savefig(outdir / f"jive_pca.1.2.png", dpi=300)

# Perform the consensus PCA
u, s, vt = scipy.sparse.linalg.svds(normalized_data, k=4)
total_var = (normalized_data**2).sum().sum()

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
        ax.set_title(format_study(study))
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.flatten()[remove_unplotted]:
        ax.remove()
    A_pct_var = s[A]**2 / total_var
    B_pct_var = s[B]**2 / total_var
    fig.supxlabel(f"{labelA} PC ({A_pct_var:0.1%})")
    fig.supylabel(f"{labelB} PC ({B_pct_var:0.1%})")
    fig.tight_layout()
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
        ax.set_title(format_study(study))
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.flatten()[remove_unplotted]:
        ax.remove()
    fig.supxlabel(f"{labelA} PC")
    fig.supylabel(f"{labelB} PC")
    fig.tight_layout()
    fig.savefig(outdir / f"individual_pca.{labelA}.{labelB}.png", dpi=300)
