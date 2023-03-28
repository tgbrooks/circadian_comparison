import pandas
import pathlib

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

jtk = pandas.read_csv(snakemake.input.jtk, sep="\t")
bootejtk = pandas.read_csv(snakemake.input.bootejtk, sep="\t")

# Load JTK and calculate overlaps
jtk = pandas.read_csv(snakemake.input.all_jtk, sep="\t")
jtk['significant'] = jtk['qvalue'] < 0.05
jtk_overlaps = []
for study1, data1 in jtk.groupby("study"):
    sig1 = set(data1.query('significant').ID)
    for study2, data2 in jtk.groupby("study"):
        if study1 == study2:
            continue
        sig2 = set(data2.query('significant').ID)
        jtk_overlaps.append({
            "study1": study1,
            "study2": study2,
            "both": len(sig1.intersection(sig2)),
            "just1": len(sig1.difference(sig2)),
            "just2": len(sig2.difference(sig1)),
        })
jtk_overlaps = pandas.DataFrame(jtk_overlaps).set_index(['study1', 'study2'])
jtk_overlaps['total'] = jtk_overlaps['both'] + jtk_overlaps['just1'] + jtk_overlaps['just2']
print(jtk_overlaps.head())
# Compute the scaling
jtk_scale = max(total.max(), jtk_overlaps['total'].max())

