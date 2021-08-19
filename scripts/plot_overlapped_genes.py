import json
import pathlib
import pandas
import numpy
import scipy.cluster.hierarchy
import matplotlib
matplotlib.use("Agg")
import pylab
from studies import sample_timepoints
import util

DPI = 300

studies = snakemake.params.studies
N_studies = len(studies)

# Load JTK values
n_genes = 1000
selected_genes = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    jtk = jtk.sort_values("ADJ.P")
    selected_genes[study] = jtk.index[:n_genes]

# Gather overlaps between each pair of studies in their top 1000 genes
results = {}
for study1 in studies: 
    study_results = {}
    for study2 in studies:
        genes1 = selected_genes[study1]
        genes2 = selected_genes[study2]
        num_common = len(genes1.intersection(genes2))
        study_results[study2] = num_common
    results[study1] = study_results

results_df = pandas.DataFrame(results)
results_df.to_csv(snakemake.output.num_common_genes, sep="\t")

# Plot these pair-wise overlaps
study_clustering = scipy.cluster.hierarchy.linkage(results_df, optimal_ordering=True)
study_ordering = scipy.cluster.hierarchy.leaves_list(study_clustering)
study_order = results_df.index[study_ordering]
fig, ax = pylab.subplots(figsize=(10,10))
h = ax.imshow(results_df.loc[study_order,study_order])
ax.xaxis.tick_top()
ax.set_xticks(numpy.arange(len(study_order)))
ax.set_xticklabels(study_order, rotation=90)
ax.set_yticks(numpy.arange(len(study_order)))
ax.set_yticklabels(study_order)
fig.colorbar(h, fraction=0.03, label="Overlap (# genes)")
fig.tight_layout()
fig.savefig(snakemake.output.num_common_genes_heatmap, dpi=DPI)

# Compute the robustness scores: number of studies where found significant
p_value_cutoff = 0.05
selected_genes = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    if study in ["Lahens15", "Manella21_Liver", "Zhang14_RNAseq_Liver_M"]:
        continue
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    #jtk = jtk.sort_values("ADJ.P")
    selected_genes[study] = jtk.index[jtk["ADJ.P"]<p_value_cutoff]

robustness_score = pandas.Series(0, index=jtk.index)

# Find intersection of all studies
intersection = selected_genes[studies[2]]
for study in selected_genes.keys():
    genes = selected_genes[study]
    intersection = intersection.intersection(genes)
    robustness_score[genes] += 1

gene_name_from_id = pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol']
intersect_DF = pandas.DataFrame({"GeneID":intersection, "GeneSymbol":intersection.map(gene_name_from_id)})
intersect_DF.to_csv(snakemake.output.common_genes_pvalue, sep="\t", index=False)
robustness_score.to_csv(snakemake.output.robustness_score, sep="\t")

# Heatmap of the genes common to all studies
data = {}
metadata = []
for study, tpmfile in zip(studies, snakemake.input.tpm):
    if study == "Weger18":
        continue
    tpm = pandas.read_csv(tpmfile, sep="\t", index_col=0)
    time = sample_timepoints(study)
    metadata.append(pandas.DataFrame({
        "study": [study for i in range(len(tpm.columns))],
        "time": time
    }, index=tpm.columns))
    # Sort TPM value columns by time
    #time_order = numpy.argsort(time)
    #tpm = tpm.iloc[:, time_order]
    data[study] = tpm

data_df = pandas.concat(data.values(), axis=1)
metadata_df = pandas.concat(metadata, axis=0)
time_order = numpy.argsort(metadata_df.time%24)
data_df = data_df.iloc[:, time_order]
heatmap = numpy.log(data_df.loc[intersection]+1)
heatmap = heatmap.divide(heatmap.median(axis=1), axis="index")
fig, ax = pylab.subplots(figsize=(60,6))
ax.imshow(heatmap, aspect = 5, interpolation='nearest')
count_study = metadata_df.groupby("study", sort=False).size().cumsum()
print(count_study)
#ax.set_xticks(count_study)
#ax.set_xticklabels(metadata_df.study.unique())
ax.set_yticks(list(range(len(intersection))))
ax.set_yticklabels(intersection)
fig.tight_layout()
fig.savefig(snakemake.output.heatmap, dpi=300)
