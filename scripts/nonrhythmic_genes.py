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

studies = snakemake.params.studies
N_studies = len(studies)

# Compute the robustness scores: number of studies where found significant
p_value_cutoff = 0.05
selected_genes = {}
jtk_amp = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    selected_genes[study] = jtk.index[(jtk["ADJ.P"]>p_value_cutoff)]
    jtk_amp[study] = jtk['AMP']

# Find intersection of all studies
intersection = selected_genes[studies[2]]
for study in selected_genes.keys():
    genes = selected_genes[study]
    intersection = intersection.intersection(genes)

#tpm mean > 0.1 check
tpm_data = {}
tpm_list = []
for study in studies:
    tpm_data[study] = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
for gene in intersection:
    flag = True
    for study in studies:
        mean = numpy.mean(tpm_data[study].loc[gene])
        if mean < 0.1: 
            flag = False
    if flag: 
        tpm_list.append(gene)

#amp/tpm mean < 0.5 check
jtk_tpm_list = []
for gene in tpm_list: 
    flag = True
    for study in studies:
        mean = numpy.mean(tpm_data[study].loc[gene])
        amp = jtk_amp[study][gene]
        value = amp/mean
        if value > 0.5: 
            flag = False
    if flag: 
        jtk_tpm_list.append(gene)

#intersect_DF = pandas.DataFrame({"GeneID":intersection, "GeneSymbol":intersection.map(gene_name_from_id)})
#intersect_DF.to_csv(snakemake.output.common_genes_pvalue, sep="\t", index=False)
pandas.Series(jtk_tpm_list).to_csv(snakemake.output.nonrhythmic_genes, sep="\t", header=None, index=None)

