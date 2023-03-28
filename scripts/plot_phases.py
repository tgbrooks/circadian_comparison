import pathlib
import numpy
import pandas
import pylab


outdir = pathlib.Path("results/Liver/phase_plots")
outdir.mkdir(exist_ok=True)

jtk = pandas.read_csv('results/Liver/jtk24.results.txt', sep='\t')
tpm = pandas.read_csv('results/Liver/tpm_all_samples.txt', sep='\t', index_col=0).drop(columns=['Symbol'])
sample_info = pandas.read_csv('results/Liver/all_samples_info.txt', sep='\t', index_col=0)
outlier_samples = [x.strip() for x in open('results/Liver/outlier_samples.txt').readlines() if x.strip()]
tpm = tpm.drop(columns=outlier_samples)

means_by_study = tpm.groupby(sample_info.study, axis=1).mean().reset_index().melt(id_vars=['Name'],value_name='mean_tpm', var_name = 'study')

data = pandas.merge(
    jtk,
    means_by_study.reset_index(),
    left_on=["ID", "study"],
    right_on=["Name", "study"],
    how="left",
).drop(columns=["Name"])
data['rel_amp'] = data.AMP / data.mean_tpm
data['rhythmic'] = data.qvalue < 0.05

high_expressed_genes = tpm.index[tpm.mean(axis=1) > 1]
sometimes_rhythmic = data.groupby("ID").rhythmic.sum() > 10
sometimes_rhythmic_genes = sometimes_rhythmic.index[sometimes_rhythmic]
N_GENES = 100
candidate_genes = set(high_expressed_genes).intersection(sometimes_rhythmic_genes)
numpy.random.seed(0)
selected_genes = numpy.random.choice(list(candidate_genes), N_GENES)

fig, ax = pylab.subplots(figsize=(12,6))
for i, gene in enumerate(selected_genes):
    gene_data = data[data.ID == gene]
    ax.scatter(
        numpy.ones(len(gene_data))*i,
        gene_data['LAG'] + numpy.random.normal(0.1, size=len(gene_data)),
        s = numpy.minimum(-numpy.log10(gene_data.qvalue), 4),
    )
ax.set_ylabel("Phase")
ax.set_xlabel("Gene")
fig.savefig(outdir / "selected_genes.phases.png", dpi=300)
