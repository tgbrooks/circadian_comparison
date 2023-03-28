import pandas
import pylab
import patsy
import scipy.stats
import numpy
import pathlib
import statsmodels.api as sm
import statsmodels.formula.api as smf

def normalize(x):
    # Set means to zero
    return (x.T - x.mean(axis=1)).T
def bh_fdr(ps):
    okay = numpy.isfinite(ps)
    qs = pandas.Series(float("NaN"), index=ps.index)
    qs[okay] = sm.stats.multipletests(ps[okay], method='fdr_bh')[1]
    return qs
def fisher(ps):
    chi2 = - 2*numpy.sum(numpy.log(ps), axis=1)
    return pandas.Series(scipy.stats.chi2.sf(chi2, df=2*ps.shape[1]), ps.index)
def bonferonni(ps):
    adj_ps = numpy.min(ps, axis=1) * numpy.sum(numpy.isfinite(ps), axis=1)
    return pandas.Series(adj_ps, index=ps.index)
def min_p(ps):
    # Assuming independence, the minimum p value follows Beta(1, n) distribution under null
    return pandas.Series(scipy.stats.beta(1,ps.shape[1]).cdf(numpy.min(ps, axis=1)), index=ps.index)

tissue = "Liver"
outdir = pathlib.Path(f"results/{tissue}/correlations")
outdir.mkdir(exist_ok=True)

tpm = pandas.read_csv(f"results/{tissue}/tpm_all_samples.txt", sep="\t", index_col=[0,1])
sample_info = pandas.read_csv(f"results/{tissue}/all_samples_info.txt", sep="\t", index_col=0)
robustness = pandas.read_csv(f"results/{tissue}/robustness_score.txt", sep="\t", index_col=0)['0']

#sample_info['time_class'] = sample_info.time.apply(lambda x: 'night' if (x % 24) >= 12 else 'day').astype('category')
sample_info['time_class'] = sample_info.time.apply(lambda x: 'night' if ((x % 24) >= 18) | ((x%24) < 6) else 'day').astype('category')

daynight_counts = (sample_info.groupby('study').time_class.value_counts().unstack())
print("Studies have the following number of timepoints per day/night class:")
print(daynight_counts)

print("Starting by-timepoint analysis")
#reps_per_timepoint = sample_info.groupby(['study', 'time']).size().groupby("study").min()
reps_per_timepoint_mod24 = sample_info.groupby(['study']).time.apply(lambda x: (x%24).value_counts().min())
num_timepoints_day = sample_info[sample_info.time_class == "day"].groupby(['study']).time.apply(lambda x : (x%24).nunique())
selected_studies = (reps_per_timepoint_mod24 >= 3) & (num_timepoints_day >= 3)
selected_studies['Janich15'] = False # Weird outlier study
print(f"We select {selected_studies.sum()} studies with at least 2 replicates per time point and at least 3 timepoints per day/night")
#selected_studies = (daynight_counts >= 8).all(axis=1)
#print(f"We select {selected_studies.sum()} studies with at least 8 replicates in each of day/night")

selected_samples = sample_info.index[sample_info.study.map(selected_studies)]
outlier_samples = set(x.strip() for x in open(f"results/{tissue}/outlier_samples.txt").readlines())
selected_samples = selected_samples.drop(outlier_samples.intersection(selected_samples))
print(f"These have {len(selected_samples)} samples")

selected_tpm = tpm[selected_samples]
selected_sample_info = sample_info.loc[selected_samples]
mean_by_study = selected_tpm.groupby(selected_sample_info.study, axis=1).mean()
#selected_genes = (mean_by_study > 30.0).all(axis=1)
#selected_genes = (mean_by_study.min(axis=1)).sort_values(ascending=False).head(300).index
selected_genes = robustness[mean_by_study.min(axis=1).values > 1.0].sort_values(ascending=False).head(200).index
selected_tpm = selected_tpm.loc[selected_genes].droplevel(1)
print(f"We choose {len(selected_tpm)} genes with high expression.")

print("Computing correlations")
results_dict = {}
for time_class, class_data in selected_tpm.groupby(selected_sample_info.time_class, axis=1):
    print(f"For {time_class}")
    results_dict[time_class] = {}
    for study, study_data in class_data.groupby(class_data.columns.map(sample_info.study), axis=1):
        print(f"For {study}")
        study_sample_info = selected_sample_info.loc[study_data.columns]
        log_tpm = numpy.log(study_data + 0.01)
        norm_log_tpm = log_tpm.groupby(study_sample_info.time, axis=1).apply(normalize)
        study_corr = norm_log_tpm.T.corr()
        results_dict[time_class][study] = study_corr
print("Done computing")
results = pandas.concat({key: pandas.concat(data) for key, data in results_dict.items()})
results.index.names = ["time_class", "study", "ID"]
results.to_csv(outdir / "day_night_correlations.txt", sep="\t")

print("Testing for differences in correlations in day versus night")
corr_differences = (results[results.index.get_level_values(0) == "day"].droplevel(0)
                        - results[results.index.get_level_values(0) == "night"].droplevel(0))
corr_differences = corr_differences.unstack(level="ID")
corr_differences = corr_differences.loc[:, corr_differences.columns.get_level_values(0) < corr_differences.columns.get_level_values(1)]

ps_r = pandas.Series(scipy.stats.ttest_1samp(corr_differences, axis=0, popmean=0)[1],
                    index=corr_differences.columns)
ps_wilcox = pandas.Series([scipy.stats.wilcoxon(diffs)[1] for i,diffs in corr_differences.T.iterrows()],
            index =  corr_differences.columns)
qs_r = bh_fdr(ps_r)
qs_wilcox = bh_fdr(ps_wilcox)

day_night_results = pandas.DataFrame({
    "p": ps_r, 
    "q": qs_r,
    "p_wilcox": ps_wilcox,
    "q_wilcox": qs_wilcox,
})
#day_night_results.to_csv(outdir / "day_night_results.txt", sep="\t")
print(f"Identified {(qs_r <0.3).sum()/2} pairs significant at q < 0.3")
est_nonnulls = len(day_night_results) -  2 * day_night_results.p.sum()
print(f"Estimate of {int(est_nonnulls)} ({est_nonnulls/len(day_night_results):0.2%}) non-null")

# By timepoint in large studies

print("Starting by-timepoint analysis")
reps_per_timepoint = sample_info.groupby(['study', 'time']).size().groupby("study").min()
#reps_per_timepoint_mod24 = sample_info.groupby(['study']).time.apply(lambda x: (x%24).value_counts().min())
selected_studies = reps_per_timepoint >= 4
print(f"Identified {sum(selected_studies)} studies with at  least 4 replicates per timepoint")

selected_samples = sample_info.index[sample_info.study.map(selected_studies)]
outlier_samples = set(x.strip() for x in open(f"results/{tissue}/outlier_samples.txt").readlines())
selected_samples = selected_samples.drop(outlier_samples.intersection(selected_samples))
print(f"These have {len(selected_samples)} samples")

selected_tpm = tpm[selected_samples]
selected_sample_info = sample_info.loc[selected_samples]
mean_by_study = selected_tpm.groupby(selected_sample_info.study, axis=1).mean()
#selected_genes = (mean_by_study.min(axis=1)).sort_values(ascending=False).head(2000).index
selected_genes = robustness.sort_values(ascending=False).head(2000).index
selected_tpm = selected_tpm.loc[selected_genes].droplevel(1)
print(f"We choose {len(selected_tpm)} genes with high expression.")

print("Beginning cosinor correlations")
corr_list = []
for study, study_data in selected_tpm.groupby(selected_sample_info.study, axis=1):
    study_sample_info = selected_sample_info.loc[study_data.columns]
    by_timepoint_corr = pandas.concat({
        time: time_data.T.corr()
            for time, time_data in study_data.groupby(study_sample_info.time, axis=1)
    }, axis=0)
    by_timepoint_corr.index.names = ["time", "ID"]
    for geneA, gene_data in by_timepoint_corr.groupby("ID"):
        time = (gene_data.index.get_level_values(0))
        corr = pandas.DataFrame({
            "r": float("NaN"),
            "cos": numpy.cos(time / 24 * 2 * numpy.pi),
            "sin": numpy.sin(time / 24 * 2 * numpy.pi),
        })
        dmatrix = patsy.dmatrix("~ cos + sin", corr)
        contrast_matrix = numpy.array([[0,1,0], [0,0,1]]) # cos=0, sin=0
        # Compute the least squares fit
        X, *_ = numpy.linalg.lstsq(dmatrix, gene_data, rcond=None)
        Y_hat = dmatrix @ X
        Y_bar = numpy.mean(gene_data, axis=0)
        N = len(corr)
        MSS = numpy.linalg.norm(Y_hat - Y_bar.values[None,:], axis=0)**2
        RSS = numpy.linalg.norm(Y_hat - gene_data, axis=0)**2
        F = (MSS/2) / (RSS/(N-3))
        p = pandas.Series(scipy.stats.f.sf(F, 2, N-3), index=gene_data.columns)
        corr_list.append(pandas.DataFrame({
            "study": study,
            "geneA": geneA,
            "geneB": p.index,
            "p": p.values,
            "amp": numpy.sqrt(X[1]**2 + X[2]**2),
            "phase": numpy.arctan2(X[2], X[1]),
            "mesor": X[0],
        }))
    print(f"Done with {study}")
corr_results = pandas.concat(corr_list, axis=0)
corr_results = corr_results[corr_results.geneA < corr_results.geneB]
corr_results.to_csv(outdir / "cosinor_correlations.txt", sep="\t", index=None)

ps = corr_results.pivot(index="study", columns=["geneA", "geneB"], values="p").astype(float)
combined_ps = fisher(ps.T)
#combined_ps = bonferonni(ps.T)
#combined_ps = min_p(ps.T)
combined_qs = bh_fdr(combined_ps)
combined_results = pandas.DataFrame({
    "p": combined_ps,
    "q": combined_qs,
})
combined_results.to_csv(outdir / "cosinor_results.txt", sep="\t")

print(f"Obtained {(combined_results.q < 0.2).sum()} pairs with q < 0.2")


# Mixed effects:
correlations_list = {}
for study, study_data in selected_tpm.groupby(selected_sample_info.study, axis=1):
    study_sample_info = selected_sample_info.loc[study_data.columns]
    by_timepoint_corr = pandas.concat({
        time: time_data.T.corr()
            for time, time_data in study_data.groupby(study_sample_info.time, axis=1)
    }, axis=0)
    by_timepoint_corr.index.names = ["time", "ID"]
    correlations_list[study] = by_timepoint_corr
correlations = pandas.concat(correlations_list, names=["study"]).reset_index().set_index("ID")
mixedlm_results_list = []
for geneA in selected_genes[:100]:
    for geneB in selected_genes[:100]:
        if geneA >= geneB:
            continue
        corr_data = correlations.loc[geneA, ['time', 'study', geneB]]
        corr_data['cos'] = numpy.cos(corr_data.time * 2 * numpy.pi / 24)
        corr_data['sin'] = numpy.sin(corr_data.time * 2 * numpy.pi / 24)
        fit = smf.mixedlm(f"{geneB} ~ cos + sin", data=corr_data, groups= "study").fit()
        mixedlm_results_list.append({
            "geneA": geneA,
            "geneB": geneB,
            'cos': fit.params['cos'],
            'sin': fit.params['sin'],
            'p_cos': fit.pvalues['cos'],
            'p_sin': fit.pvalues['sin'],
            'converge': fit.converged,
        })
mixedlm_results = pandas.DataFrame(mixedlm_results_list)

def plot_corr(geneA, geneB):
    #fig, axes = pylab.subplots(nrows=selected_sample_info.study.nunique())
    fig, ax = pylab.subplots()
    for (study, study_data) in (selected_tpm.groupby(selected_sample_info.study, axis=1)):
        study_sample_info = selected_sample_info.loc[study_data.columns]
        by_timepoint_corr = pandas.Series({
            time: time_data.T[[geneA, geneB]].corr().iloc[0,1]
                for time, time_data in study_data.groupby(study_sample_info.time, axis=1)
        })
        by_timepoint_corr.index.names = ["time"]

        #ax.scatter(by_timepoint_corr.index%24, by_timepoint_corr.values)
        #ax.set_title(study)
        ax.plot(by_timepoint_corr.index, by_timepoint_corr.values, label=study)
        #ax.set_xlim(0,24)
        ax.set_xticks(numpy.arange(0,36,6))
        ax.set_xlabel("ZT")
        ax.set_ylim(-1,1)

def plot_daynight(geneA, geneB):
    fig, axes = pylab.subplots(figsize=(16,4), ncols=selected_studies.sum(), nrows=2)
    for row, (time_class, class_data) in zip(axes, selected_tpm.groupby(selected_sample_info.time_class, axis=1)):
        for ax, (study, study_data) in zip(row, class_data.groupby(class_data.columns.map(sample_info.study), axis=1)):
            study_sample_info = selected_sample_info.loc[study_data.columns]
            log_tpm = numpy.log(study_data + 0.01)
            norm_log_tpm = log_tpm.groupby(study_sample_info.time, axis=1).apply(normalize)
            ax.scatter(norm_log_tpm.loc[geneA], norm_log_tpm.loc[geneB], label=study)
            ax.set_aspect(1)
            ax.set_xticks([])
            ax.set_yticks([])
            if ax in axes[0]:
                ax.set_title(study)

    fig.tight_layout()


