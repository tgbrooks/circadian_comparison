import pandas
import statsmodels.stats

from studies import studies

jtk = {}
mean_reads = {}
for study in studies:
    this_tpm = pandas.read_csv(f"data/{study}/expression.num_reads.txt", sep="\t", index_col=0)
    this_jtk = pandas.read_csv(f"data/{study}/jtk/JTKresult_expression.tpm.txt", sep="\t", index_col=0)
    this_jtk = this_jtk.loc[this_tpm.index]
    jtk[study] = this_jtk
    mean_reads[study] = this_tpm.mean(axis=1)

Q_CUTOFF = 0.10
mean_reads_cutoffs = [0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25]

num_significant = {}
for study in studies:
    num_significant[study] = {}
    for cutoff in mean_reads_cutoffs:
        selected = mean_reads[study] >= cutoff
        ps = jtk[study].loc[selected, 'ADJ.P']
        _, qs, _, _ = statsmodels.stats.multitest.multipletests(ps, method="fdr_bh")

        num_significant[study][cutoff] = (qs < Q_CUTOFF).sum()

num_significant = pandas.DataFrame(num_significant)
num_significant.to_csv("results/jtk.num_signifcant.by_cutoffs.txt", sep="\t")
