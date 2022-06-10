import pathlib
import pandas

from studies import targets

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

for count_type, count_file in zip(['tpm', 'num_reads'], [snakemake.input.tpm, snakemake.input.num_reads]):
        count_data = pandas.read_csv(count_file, sep="\t")
        count_data.to_csv(outdir / f"{count_type}.by_sample.txt.gz", sep="\t")

(outdir / "study_metadata.txt").write_text(pathlib.Path(snakemake.input.study_table).read_text())

outlier_samples = [x.strip() for x in open(snakemake.input.outliers).readlines()]

sample_metadata = pandas.read_csv(snakemake.input.sample_info, sep="\t")
sample_metadata.columns = ['sample', 'study', 'time']
sample_metadata['outlier'] = sample_metadata['sample'].isin(outlier_samples)
sample_metadata['study'] =  sample_metadata.study.map(lambda x: targets[x]['short_name'])
sample_metadata.to_csv(outdir / "sample_metadata.txt", index=None, sep="\t")

jtk = pandas.read_csv(snakemake.input.jtk, sep="\t")
jtk['study'] =  jtk.study.map(lambda x: targets[x]['short_name'])
jtk.to_csv(outdir / "jtk.results.txt.gz", index=None, sep="\t")

robustness_score = pandas.read_csv(snakemake.input.robustness_score, sep="\t")
robustness_score.columns = ['ID', 'score']
robustness_score.to_csv(outdir / "robustness_score.txt", index=None, sep="\t")

simdir = outdir / "SIM"
simdir.mkdir(exist_ok=True)

(simdir / "summary.txt").write_text(pathlib.Path(snakemake.input.sim_summary).read_text())
(simdir / "curves_fit.txt").write_text(pathlib.Path(snakemake.input.sim_curves_fit).read_text())
(simdir / "curves_pstd.txt").write_text(pathlib.Path(snakemake.input.sim_curves_pstd).read_text())
(simdir / "re_structure.txt").write_text(pathlib.Path(snakemake.input.sim_re_structure).read_text())
sim_re = pandas.read_csv(snakemake.input.sim_re, sep="\t")
sim_re['study'] =  sim_re.study.map(lambda x: targets[x]['short_name'])
sim_re.to_csv(simdir / "re.txt", sep="\t", index=None)
