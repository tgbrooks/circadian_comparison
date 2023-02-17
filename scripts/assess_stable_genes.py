import pandas

# Load basic information
tpm = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=0).drop(columns=["Symbol"])
sample_info = pandas.read_csv(snakemake.input.sample_info, sep="\t", index_col=0)
outlier_samples = [x.strip() for x in open(snakemake.input.outliers).readlines()]
sample_info = sample_info[~sample_info.index.isin(outlier_samples)]

# Summarize TPM
mean_tpm_by_study = tpm.groupby(sample_info.study, axis=1).mean()
std_tpm_by_study = tpm.groupby(sample_info.study, axis=1).std()
rel_std_by_study = std_tpm_by_study / mean_tpm_by_study

# Load spline fits
spline_fit_summary = pandas.read_csv(snakemake.input.spline_fit_summary, sep="\t", index_col=0)

# Load JTK values
jtk = pandas.read_csv(snakemake.input.jtk, sep="\t", index_col=0)

# Determine stable genes
info = pandas.DataFrame({
    "min_mean_tpm": mean_tpm_by_study.min(axis=1),
    "max_rel_std": rel_std_by_study.max(axis=1),
    "min_jtk_q": jtk.groupby("ID").qvalue.min(),
    "spline_fit_rhythmic": spline_fit_summary.is_rhythmic,
})
print(info.head())
info['stable'] = (info.min_mean_tpm > 1) & (info.max_rel_std < 0.5) & (info.min_jtk_q > 0.05) & (info.spline_fit_rhythmic != True)

print(f"Identified {info.stable.sum()} stable genes out of {len(info)}")
print(f"20 sampled stable genes:\n{info.query('stable').sample(20)}")

pandas.Series(info.query('stable').index).to_csv(snakemake.output.stable_gene_list, sep="\t", index=None, header=None)
info.sort_index().to_csv(snakemake.output.stable_table, sep="\t")
