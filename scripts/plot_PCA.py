import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab
import statsmodels.api as sm
from studies import sample_timepoints

from styles import color_by_study, shape_by_study

studies = snakemake.params.studies
DPI = 300
PSEUDOCOUNT = 0.01

# Ensure our output directory exists
OUT_DIR = pathlib.Path(snakemake.output.dir)
OUT_DIR.mkdir(exist_ok=True)

data = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=[0,1])
sample_info = pandas.read_csv(snakemake.input.sample_info, sep="\t", index_col=[0])
outlier_samples = [x.strip() for x in open(snakemake.input.outlier_samples).readlines()]

jtk_amp = {}
jtk_p = {}
for study, jtk_file in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtk_file, sep="\t", index_col=0).sort_index()
    #jtk.loc[jtk.dropped, 'ADJ.P'] = 1 # JTK gives 0 p-values to bad data: we'll use 1 instead
    jtk_amp[study] = jtk.AMP
    jtk_p[study] = jtk['ADJ.P']

with open("results/study_classification.json") as f:
    study_classification = json.load(f)

robustness = pandas.read_csv(snakemake.input.robustness, sep="\t", index_col=0)['0'].sort_values(ascending=False)

jtk_amp_df = pandas.DataFrame(jtk_amp)
jtk_p_df = pandas.DataFrame(jtk_p)

# Drop all rows with no non-zero values
data = data[~(data == 0).all(axis=1)]
jtk_amp_df = jtk_amp_df[~(jtk_amp_df == 0).all(axis=1)]
jtk_p_df = jtk_p_df[~(jtk_p_df == 1).all(axis=1)]

#Drop outlier samples
data = data.drop(columns=outlier_samples)
sample_info = sample_info.drop(index=outlier_samples)
print(f"data has shape {data.shape}")
print(f"sample_info has shape {sample_info.shape}")
print(f"jtk_amp has shape {jtk_amp_df.shape}")

#Compute PCA on the expression levels (TPM)
used_studies = sample_info.study.unique()
log_data = numpy.log(data+PSEUDOCOUNT)
PCA_by_TPM= sm.PCA(log_data.T, ncomp=2, method="nipals")

# Color by study
fig, ax = pylab.subplots(figsize=(12,12))
for study, metadata in sample_info.groupby('study'):
    ax.scatter(
        PCA_by_TPM.factors.loc[metadata.index, "comp_0"],
        PCA_by_TPM.factors.loc[metadata.index, "comp_1"],
        color=color_by_study[study],
        marker=shape_by_study[study],
        label=study,
    )
ax.set_xlabel(f"PCA 1 ({PCA_by_TPM.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA_by_TPM.rsquare[2] - PCA_by_TPM.rsquare[1]:0.0%})")
fig.legend(fontsize = 'x-small', ncol=2)
fig.savefig(OUT_DIR/f"all_samples_study.TPM.png", dpi=DPI)
fig.savefig(OUT_DIR/f"all_samples_study.TPM.svg")


# Color by time-of-day
colors=sample_info.time%24

fig, ax = pylab.subplots(figsize=(12,12))
h=ax.scatter(
        PCA_by_TPM.factors["comp_0"],
        PCA_by_TPM.factors["comp_1"],
        c=colors,
        )
ax.set_xlabel(f"PCA 1 ({PCA_by_TPM.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA_by_TPM.rsquare[2] - PCA_by_TPM.rsquare[1]:0.0%})")
fig.colorbar(h)
fig.savefig(OUT_DIR/f"all_samples_time.TPM.png", dpi=DPI)
fig.savefig(OUT_DIR/f"all_samples_time.TPM.svg")


# Color by study pca_type (paired/stranded/etc)
fig, ax = pylab.subplots(figsize=(12,12))
sample_info['study_seq_type'] = sample_info.study.map(study_classification)
for study, metadata in sample_info.groupby('study_seq_type'):
    ax.scatter(
        PCA_by_TPM.factors.loc[metadata.index, "comp_0"],
        PCA_by_TPM.factors.loc[metadata.index, "comp_1"],
        label=study,
    )
ax.set_xlabel(f"PCA 1 ({PCA_by_TPM.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA_by_TPM.rsquare[2] - PCA_by_TPM.rsquare[1]:0.0%})")
fig.legend(ncol=2)
fig.savefig(OUT_DIR/f"all_samples_study_classification.TPM.png", dpi=DPI)
fig.savefig(OUT_DIR/f"all_samples_study_classification.TPM.svg")

## Compute PCA by the JTK results (amplitude)
# One point will be one study, not one sample
log_amp_data = numpy.log(jtk_amp_df+PSEUDOCOUNT)
log_amp_data = log_amp_data[log_amp_data.index.isin(robustness[:1000].index)]
PCA_by_JTK_amp = sm.PCA(log_amp_data.T, ncomp=2, method="nipals")

# Color by study
fig, ax = pylab.subplots(figsize=(5,5))
for study, amp in log_amp_data.iteritems():
    ax.scatter(
        PCA_by_JTK_amp.factors.loc[study, "comp_0"],
        PCA_by_JTK_amp.factors.loc[study, "comp_1"],
        color=color_by_study[study],
        marker=shape_by_study[study],
        label=study,
    )
ax.set_xlabel(f"PCA 1 ({PCA_by_JTK_amp.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA_by_JTK_amp.rsquare[2] - PCA_by_JTK_amp.rsquare[1]:0.0%})")
fig.legend(fontsize = 'x-small', ncol=2)
fig.savefig(OUT_DIR/f"all_studies.JTK_amp.png", dpi=DPI)
fig.savefig(OUT_DIR/f"all_studies.JTK_amp.svg")

# Color by sequencing classification
fig, ax = pylab.subplots(figsize=(5,5))
classifications = sorted(str(x) for x in sample_info.study_seq_type.unique())
for classification in classifications:
    studies = [study for study, class_ in study_classification.items()
                    if class_ == classification
                        and study in jtk_p_df.columns]
    ax.scatter(
        PCA_by_JTK_amp.factors.loc[studies, "comp_0"],
        PCA_by_JTK_amp.factors.loc[studies, "comp_1"],
        label=classification,
    )
ax.set_xlabel(f"PCA 1 ({PCA_by_JTK_amp.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA_by_JTK_amp.rsquare[2] - PCA_by_JTK_amp.rsquare[1]:0.0%})")
fig.legend(ncol=2)
fig.savefig(OUT_DIR/f"all_studies.classification.JTK_amp.png", dpi=DPI)
fig.savefig(OUT_DIR/f"all_studies.classification.JTK_amp.svg")


## Compute PCA by the JTK results (p-value)
# One point will be one study, not one sample
log_p_data = numpy.log(jtk_p_df)
log_p_data = log_p_data[log_p_data.index.isin(robustness[:1000].index)]
PCA_by_JTK_p= sm.PCA(log_p_data.T, ncomp=2, method="nipals")

# Color by study
fig, ax = pylab.subplots(figsize=(5,5))
for study, amp in log_p_data.iteritems():
    ax.scatter(
        PCA_by_JTK_p.factors.loc[study, "comp_0"],
        PCA_by_JTK_p.factors.loc[study, "comp_1"],
        color=color_by_study[study],
        marker=shape_by_study[study],
        label=study,
    )
ax.set_xlabel(f"PCA 1 ({PCA_by_JTK_p.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA_by_JTK_p.rsquare[2] - PCA_by_JTK_p.rsquare[1]:0.0%})")
fig.legend(fontsize = 'x-small', ncol=2)
fig.savefig(OUT_DIR/f"all_studies.JTK_p.png", dpi=DPI)
fig.savefig(OUT_DIR/f"all_studies.JTK_p.svg")
