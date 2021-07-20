import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

import util

studies = snakemake.params.studies
N_studies = len(studies)

n_genes = 1000
selected_genes = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    jtk = jtk.sort_values("ADJ.P")
    selected_genes[study] = jtk.index[:n_genes]

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

p_value_cutoff = 0.05
selected_genes = {}
for study, jtkfile in zip(studies, snakemake.input.jtk):
    jtk = pandas.read_csv(jtkfile, sep="\t", index_col=0)
    #jtk = jtk.sort_values("ADJ.P")
    selected_genes[study] = jtk.index[jtk["ADJ.P"]<p_value_cutoff]

intersection = selected_genes[studies[0]]
for study in studies:
    genes = selected_genes[study]
    intersection = intersection.intersection(genes)

gene_name_from_id = pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol']
intersect_DF = pandas.DataFrame({"GeneID":intersection, "GeneSymbol":intersection.map(gene_name_from_id)})
intersect_DF.to_csv(snakemake.output.common_genes_pvalue, sep="\t", index=False)
