import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

studies = snakemake.params.studies
N_studies = len(studies)

READ_UNITS = 1e6 # Use millions of reads
DPI = 300

salmon_metainfos = [f"data/{study}/salmon.meta_info.json" for study in studies]
mapping_percents = []
num_mapped = []
num_processed = []
for metainfo_file in salmon_metainfos:
    with open(metainfo_file) as stream:
        metainfo = json.load(stream)
    mapping_percents.append([x['percent_mapped'] for k,x in metainfo.items()])
    num_mapped.append([x['num_mapped']/READ_UNITS for k,x in metainfo.items()])
    num_processed.append([x['num_processed']/READ_UNITS for k,x in metainfo.items()])

fig, ax = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
ax.boxplot(mapping_percents)
ax.set_xticks(numpy.arange(N_studies)+1)
ax.set_xticklabels(studies, rotation=90)
ax.set_ylabel("Percent Mapped")
fig.tight_layout()
fig.savefig(snakemake.output.percent_mapping, dpi=DPI)

fig, ax = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
ax.boxplot(num_mapped)
ax.set_xticks(numpy.arange(N_studies)+1)
ax.set_xticklabels(studies, rotation=90)
ax.set_ylabel("Num Reads Mapped (M)")
ax.set_yscale("log")
fig.tight_layout()
fig.savefig(snakemake.output.num_mapped, dpi=DPI)

fig, ax = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
ax.boxplot(num_processed)
ax.set_xticks(numpy.arange(N_studies)+1)
ax.set_xticklabels(studies, rotation=90)
ax.set_ylabel("Num Reads Total (M)")
ax.set_yscale("log")
fig.tight_layout()
fig.savefig(snakemake.output.num_processed, dpi=DPI)

# Plot certain genes for QC
genes_ID = ["ENSMUSG00000086503","ENSMUSG00000029368", "ENSMUSG00000064337", "ENSMUSG00000064339", "ENSMUSG00000064336"]
genes_symbol = ["Xist", "Alb", "mt-Rnr1", "mt-Rnr2", "mt-Tf"]
tpm = {}
for study in studies:
    tpm[study] = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
for gene, name in zip(genes_ID, genes_symbol):
    fig, ax = pylab.subplots(figsize=(0.7*N_studies+0.7,4))
    expression_data = []
    for study in studies:
        expression_data.append(tpm[study].loc[gene])
    ax.boxplot(expression_data)

    ax.set_xticks(numpy.arange(N_studies)+1)
    ax.set_xticklabels(studies, rotation=90)
    ax.set_ylabel(name + " Expression")
    ax.set_title(gene)
    fig.tight_layout()
    fig.savefig(snakemake.output[gene], dpi=DPI)

