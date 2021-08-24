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

data = {}
metadata = []
for study, tpmfile in zip(studies, snakemake.input.tpm):
    if study == "Weger18": 
        continue
    tpm = pandas.read_csv(tpmfile, sep="\t", index_col=0)
    data[study] = tpm
    metadata.append(pandas.DataFrame({
        "study": [study for i in range(len(tpm.columns))], 
        "time": sample_timepoints(study), 
        }, index=tpm.columns))

with open("results/study_classification.json") as f:
    study_classification = json.load(f)

data_df = pandas.concat(data.values(), axis=1)
metadata_df = pandas.concat(metadata, axis=0)
# Drop all rows with no non-zero values
data_df = data_df[~(data_df == 0).all(axis=1)]
print(f"data has shape {data_df.shape}")

#Compute PCA
used_studies = metadata_df.study.unique()
log_data = numpy.log(data_df+1)
PCA = sm.PCA(log_data.T, ncomp=2, method="nipals")

# Color by study
fig, ax = pylab.subplots(figsize=(12,12))
for study, metadata in metadata_df.groupby('study'):
    ax.scatter(
        PCA.factors.loc[metadata.index, "comp_0"],
        PCA.factors.loc[metadata.index, "comp_1"],
        color=color_by_study[study],
        marker=shape_by_study[study],
        label=study,
    )
ax.set_xlabel(f"PCA 1 ({PCA.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA.rsquare[2] - PCA.rsquare[1]:0.0%})")
fig.legend(fontsize = 'x-small', ncol=2)
fig.savefig(snakemake.output.all_samples_study, dpi=DPI)
fig.savefig(snakemake.output.all_samples_study_svg)


# Color by time-of-day
colors=metadata_df.time%24

fig, ax = pylab.subplots(figsize=(12,12))
h=ax.scatter(
        PCA.factors["comp_0"],
        PCA.factors["comp_1"],
        c=colors,
        )
ax.set_xlabel(f"PCA 1 ({PCA.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA.rsquare[2] - PCA.rsquare[1]:0.0%})")
fig.colorbar(h)
fig.savefig(snakemake.output.all_samples_time, dpi=DPI)
fig.savefig(snakemake.output.all_samples_time_svg)


# Color by study type (paired/stranded/etc)
fig, ax = pylab.subplots(figsize=(12,12))
metadata_df['study_type'] = metadata_df.study.map(study_classification)
for study, metadata in metadata_df.groupby('study_type'):
    ax.scatter(
        PCA.factors.loc[metadata.index, "comp_0"],
        PCA.factors.loc[metadata.index, "comp_1"],
        label=study,
    )
ax.set_xlabel(f"PCA 1 ({PCA.rsquare[1]:0.0%})")
ax.set_ylabel(f"PCA 2 ({PCA.rsquare[2] - PCA.rsquare[1]:0.0%})")
fig.legend(fontsize = 'x-small', ncol=2)
fig.savefig(snakemake.output.all_samples_study_classification, dpi=DPI)
fig.savefig(snakemake.output.all_samples_study_classification_svg)

