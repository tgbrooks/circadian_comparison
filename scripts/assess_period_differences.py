import pathlib

import scipy.stats
import pandas
import numpy
import pylab

import studies

tissue = snakemake.wildcards.tissue
outdir = pathlib.Path(f"results/{tissue}/assess_jtk")
outdir.mkdir(exist_ok=True)

jtk = pandas.read_csv(f"results/{tissue}/jtk.results.txt", sep="\t", index_col=[0,1])
robustness = pandas.read_csv(f"results/{tissue}/robustness_score.txt", sep="\t", index_col=0)['0'].sort_index()

# Find studies with at least 8 timepoints per day, since these have the resolution to detect
# periods other than 24. And we require all to be 'highlighted' studies, i.e. from consistent
# conditions and so comparable
sample_info = pandas.read_csv(f"results/{tissue}/all_samples_info.txt", sep="\t", index_col=0)
num_timepoints = sample_info.groupby("study").time.apply(lambda x: (x%24).nunique())
highres_study = num_timepoints.index[num_timepoints >= 8]
selected_study = [study for study in highres_study if studies.targets[study].get("highlight", False)]
if len(selected_study) < 2:
    print("No suitable studies. Aborting")
    (outdir / "period_statistics.txt").touch()
    exit()

# Select just data when JTK significant and from high resolution datasets
sig = (~jtk.dropped) & (jtk['qvalue'] < 0.05)
good_jtk = jtk[sig]
period = good_jtk.reset_index().pivot(index="ID", columns="study", values="PER")[highres_study]
period = period.sort_index()

# Compute the Cramer's V and Spearman correlation statistics between the selected studies
statistics_list = []
for studyA in selected_study:
    for studyB in selected_study:
        if studyA == studyB:
            continue
        contingency = pandas.crosstab(period[studyA], period[studyB])
        chi2, p, ddof, expected = scipy.stats.chi2_contingency(contingency)
        cramerV = numpy.sqrt(chi2 / (contingency.sum().sum() * (numpy.min(contingency.shape) - 1)))
        spearmanr, p = scipy.stats.spearmanr(period[studyA], period[studyB], nan_policy="omit")
        statistics_list.append({
            "studyA": studyA,
            "studyB": studyB,
            "CramerV": cramerV,
            "SpearmanR": spearmanr,
        })
statistics = pandas.DataFrame(statistics_list)
statistics.to_csv(outdir / "period_statistics.txt", sep="\t", index=False)
print("Median PER statistics comparing pairs of {len(selected_study)} studies:")
print(statistics.median(axis=0))
print("Max PER statistics comparing pairs of {len(selected_study)} studies:")
print(statistics[['CramerV', 'SpearmanR']].max(axis=0))

# Check if any genes are consistently non-24
# Just look at < 24 since few studies can detect >24
# Require no >=24 periods and at least four <24 hours periods among those significant
def find_low_period(period):
    return period[ ((period < 24) | (period.isna())).all(axis=1) & (period.count(axis=1) >= 4)]
low_period = find_low_period(period)
print(f"Identified {len(low_period)} low-period genes")
low_period.to_csv(outdir/ "low_period_genes.txt", sep="\t")
random_num_low_period = []
def permute_finite(x):
    ''' permutation of finite elements of x, leaving nans unchanged '''
    out = x.copy()
    finite = ~pandas.isna(x)
    out[finite] = numpy.random.permutation(out[finite])
    return out
for i in range(100):
    shuffled_period = pandas.DataFrame({key: permute_finite(col.values)
                                        for key,col in period.items()})
    random_num_low_period.append(len(find_low_period(shuffled_period)))
print("Random permutations give low periods:")
print(pandas.Series(random_num_low_period).describe())

## By robustness
# Categories of robustness
high = (robustness >= 35)
mid = (robustness < 35) & (robustness >= 20)
low = (robustness < 20)

# Period histograms by robustness
bins = numpy.arange(19,30)
fig, ax = pylab.subplots(figsize=(4,4))
ax.hist(period[robustness >= 35].unstack(), bins=bins, density=True, alpha=0.5, label=">=35")
ax.hist(period[(robustness < 35) & (robustness >= 20)].unstack(), bins=bins, density=True, alpha=0.5, label="20-34")
ax.hist(period[(robustness < 20)].unstack(), bins=bins, density=True, alpha=0.5, label="0-19")
ax.set_xlabel("JTK Period")
ax.legend()
fig.savefig(outdir/"period_distributions.png", dpi=300)

# Tabular form: period by robustness
counts_by_period = pandas.DataFrame({
    ">=35": period[high].unstack().value_counts() / period[high].count().sum(),
    "20-34": period[mid].unstack().value_counts() / period[mid].count().sum(),
    "0-19": period[low].unstack().value_counts() / period[low].count().sum(),
}).fillna(0)
counts_by_period.index.name = "period"
counts_by_period.index = counts_by_period.index.astype(int)
(counts_by_period*100).to_csv(outdir / "period_counts.txt", sep="\t", float_format="%0.1f")

# Phase differences
#
## Compute the differences between studies
## Using two ways:
## 1. mean absolute difference of phase (done cyclicly, so 24=0)
## 2. number of highly discordant genes (per difference > 6 hrs)
#mad_dict = {}
#discordant_dict = {}
#for study1, jtk1 in jtk.groupby("study"):
#    mad_dict[study1] = {}
#    discordant_dict[study1] = {}
#    for study2, jtk2 in jtk.groupby("study"):
#        in_both = (jtk1.qvalue < 0.1) & (jtk2.qvalue < 0.1)
#        #mad = scipy.stats.circmean(jtk1[in_both].LAG - jtk2[in_both].LAG, 0, 24)
#        diff = (jtk1.LAG[in_both] - jtk2.LAG[in_both])%24
#        diff[diff > 12] = 24 - diff[diff > 12]
#        mad_dict[study1][study2]  = diff.mean()
#        discordant_dict[study1][study2] = sum(diff >= 6)
#mean_abs_diff = pandas.DataFrame(mad_dict)
#num_discordant_by_study = pandas.DataFrame(discordant_dict)
#
#
#def circ_abs_diff(x,y):
#    ''' circular absolute difference between two phases '''
#    diff = (x - y) % 24
#    if diff > 12:
#        return 24 - diff
#    else:
#        return diff
#
#def count_discordant_phases(phases):
#    ''' number of phase pairs that are discordant (>= 6 hrs apart) '''
#    return len([x for x in phases for y in phases if circ_abs_diff(x,y) >= 6 and x > y])
#
#num_discordant_by_gene = jtk.groupby(jtk.index).apply(lambda x: count_discordant_phases(x[x.qvalue < 0.1].LAG))
