import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab
import statsmodels.api as sm
import util
from studies import sample_timepoints

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

data_df = pandas.concat(data.values(), axis=1)
metadata_df = pandas.concat(metadata, axis=0)
# Drop all rows with no non-zero values
data_df = data_df[~(data_df == 0).all(axis=1)]
print(f"data has shape {data_df.shape}")

color_by_study = dict(zip(metadata_df.study.unique(), [pylab.get_cmap("tab10")(i) for i in range(10)]))
colors = metadata_df.study.map(color_by_study)
log_data = numpy.log(data_df+1)
PCA = sm.PCA(log_data.T, ncomp=2, method="nipals")

fig, ax = pylab.subplots(figsize=(12,12))
ax.scatter(PCA.factors["comp_0"], PCA.factors["comp_1"], c=colors)
util.legend_from_colormap(fig, color_by_study)
fig.savefig(snakemake.output.all_samples_study, dpi=DPI)

colors=metadata_df.time%24

fig, ax = pylab.subplots(figsize=(12,12))
h=ax.scatter(PCA.factors["comp_0"], PCA.factors["comp_1"], c=colors)
fig.colorbar(h)
fig.savefig(snakemake.output.all_samples_time, dpi=DPI)
