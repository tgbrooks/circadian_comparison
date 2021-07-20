import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab
import statsmodels.api as sm
import util

studies = snakemake.params.studies
DPI = 300

data = {}
for study, tpmfile in zip(studies, snakemake.input.tpm):
    tpm = pandas.read_csv(tpmfile, sep="\t", index_col=0)
    data[study] = tpm

data_df = pandas.concat(data.values(), axis=1)
# Drop all rows with no non-zero values
data_df = data_df[~(data_df == 0).all(axis=1)]
print(f"data has shape {data_df.shape}")

log_data = numpy.log(data_df+1)
PCA = sm.PCA(log_data.T, ncomp=2, method="nipals")

fig, ax = pylab.subplots(figsize=(8,8))
ax.scatter(PCA.factors["comp_0"], PCA.factors["comp_1"])
fig.savefig(snakemake.output.all_samples, dpi=DPI)

