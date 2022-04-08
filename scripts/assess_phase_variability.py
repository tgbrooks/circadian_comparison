import pandas
import pylab
import seaborn as sns
import scipy.special

PHASE_STD_CUTOFF = 2 # in hours

summary = pandas.read_csv(snakemake.input.spline_fit_summary, sep="\t", index_col=0)

# Convert phi values into phase
# NOTE: this isn't realy a standard deviation of phase values due to the fact
# that this transformation is non-linear and so f(std(x)) != std(f(x)) in general
summary['phase_std'] = scipy.special.expit(summary.re_phi_sd) * 24 - 12

fig, ax= pylab.subplots(figsize=(3,3), constrained_layout=True)
sns.histplot(
    x = "phase_std",
    data = summary.query("is_rhythmic"),
    bins = 12,
    ax = ax,
)
ax.set_xlabel("Phase Variability (hours)")
ax.set_ylabel("Num. Genes")
fig.savefig(snakemake.output.phase_std_distribution, dpi=300)

print(f"Identified {(summary.query('is_rhythmic').phase_std > PHASE_STD_CUTOFF).sum()} with phase STD above {PHASE_STD_CUTOFF}")
print(f"Corresponding re_phi_std is {scipy.special.logit((PHASE_STD_CUTOFF + 12)/24)}")
print(f"Mean value {summary.query('is_rhythmic').phase_std.mean()}")
