import json
import pathlib
import pandas
import numpy
import scipy.cluster.hierarchy
import matplotlib
matplotlib.use("Agg")
import pylab
from studies import sample_timepoints, targets
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
    selected_genes[study] = jtk[~jtk.dropped].index[:n_genes] #Top genes by JTK

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

# Overlaps by q-value / p-value
jtk = pandas.read_csv(snakemake.input.jtk_results, sep="\t")
jtk['q_significant'] = jtk.qvalue < 0.05
jtk['p_significant'] = jtk['ADJ.P'] < 0.05
qvalues_overlap = []
for studyA, dataA in jtk.groupby("study"):
    for studyB, dataB in jtk.groupby("study"):
        joined = pandas.merge(
            dataA,
            dataB,
            left_on="ID",
            right_on="ID",
            suffixes=("_A", "_B")
        )
        qvalues_overlap.append({
            "A": studyA,
            "B": studyB,
            "A_q_significant": dataA.q_significant.sum(),
            "A_q_B_p_significant": (joined.q_significant_A & joined.p_significant_B).sum(),
            "A_q_B_q_significant": (joined.q_significant_A & joined.q_significant_B).sum(),
        })
qvalues_overlap_df = pandas.DataFrame(qvalues_overlap)
qvalues_overlap_df.to_csv(snakemake.output.q_values_overlaps, sep="\t")

# Plot these pair-wise overlaps
study_clustering = scipy.cluster.hierarchy.linkage(results_df, optimal_ordering=True)
study_ordering = scipy.cluster.hierarchy.leaves_list(study_clustering)
study_order = results_df.index[study_ordering]
fig, ax = pylab.subplots(figsize=(8,8))
h = ax.imshow(results_df.loc[study_order,study_order], vmin=0, vmax=n_genes)
ax.xaxis.tick_top()
ax.set_xticks(numpy.arange(len(study_order)))
ax.set_xticklabels([targets[x]['short_name'] for x in study_order], rotation=90)
ax.set_yticks(numpy.arange(len(study_order)))
ax.set_yticklabels([targets[x]['short_name'] for x in study_order])
fig.colorbar(h, fraction=0.03, label="Overlap (# genes)")
fig.tight_layout()
fig.savefig(snakemake.output.num_common_genes_heatmap, dpi=DPI)

# Plot these pair-wise overlaps of q-values
qvalue_results = qvalues_overlap_df.pivot("A", "B", "A_q_B_q_significant")
fig, ax = pylab.subplots(figsize=(8,8))
h = ax.imshow(qvalue_results.loc[study_order,study_order], vmin=0, vmax=n_genes)
ax.xaxis.tick_top()
ax.set_xticks(numpy.arange(len(study_order)))
ax.set_xticklabels([targets[x]['short_name'] for x in study_order], rotation=90)
ax.set_yticks(numpy.arange(len(study_order)))
ax.set_yticklabels([targets[x]['short_name'] for x in study_order])
fig.colorbar(h, fraction=0.03, label="Overlap (# genes)")
fig.tight_layout()
fig.savefig(snakemake.output.num_qvalue_overlap_heatmap, dpi=DPI)

# Plot these pair-wise overlaps of q-values vs p-values
qvalue_results = qvalues_overlap_df.pivot("A", "B", "A_q_B_p_significant")
fig, ax = pylab.subplots(figsize=(8,8))
h = ax.imshow(qvalue_results.loc[study_order,study_order], vmin=0, vmax=n_genes)
ax.xaxis.tick_top()
ax.set_xticks(numpy.arange(len(study_order)))
ax.set_xticklabels([targets[x]['short_name'] for x in study_order], rotation=90)
ax.set_yticks(numpy.arange(len(study_order)))
ax.set_yticklabels([targets[x]['short_name'] for x in study_order])
fig.colorbar(h, fraction=0.03, label="Overlap (# genes)")
fig.tight_layout()
fig.savefig(snakemake.output.num_q_p_value_overlap_heatmap, dpi=DPI)

# Compute the robustness scores: number of studies where found significant
p_value_cutoff = 0.05
selected_genes = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    selected_genes[study] = jtk.index[(jtk["ADJ.P"]<p_value_cutoff) & (~jtk.dropped)]

# Find intersection and robustness scores of all studies
robustness_score = pandas.Series(0, index=jtk.index)
intersection = selected_genes[studies[2]]
for study in selected_genes.keys():
    genes = selected_genes[study]
    intersection = intersection.intersection(genes)
    robustness_score[genes] += 1
robustness_score.to_csv(snakemake.output.robustness_score, sep="\t")

gene_name_from_id = pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol']
intersect_DF = pandas.DataFrame({"GeneID":intersection, "GeneSymbol":intersection.map(gene_name_from_id)})
intersect_DF.to_csv(snakemake.output.common_genes_pvalue, sep="\t", index=False)
