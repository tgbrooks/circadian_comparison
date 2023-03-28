import json

import pandas
import numpy
import statsmodels.formula.api as smf

import studies
PSEUDOCOUNT = 0.01

tpm = pandas.read_csv("results/Liver/tpm_all_samples.txt", sep="\t", index_col=[0]).select_dtypes(['number'])

sample_info = pandas.read_csv("results/Liver/all_samples_info.txt", sep="\t", index_col=0)

with open("results/study_classification.json") as f:
    study_classification = json.load(f)
seq_type = sample_info.study.map(lambda study: studies.targets[study].get('seq', 'PolyA'))
sample_info['seq_type'] = seq_type + "-" + sample_info.study.map(study_classification)

robustness = pandas.read_csv("results/Liver/robustness_score.txt", sep="\t", index_col=0)
top_circadian = robustness.sort_values(ascending=False)[:1000]

# Now we partition the variance up by several factors:
# 1. Sequencing type
# 2. Study (random effects?)
# 3. Time of day
for gene in top_circadian:
    d = pandas.DataFrame({
        "tpm": numpy.log(tpm.loc[gene]+PSEUDOCOUNT),
        "seq_type": tpm.columns.map(sample_info.seq_type),
        "study": tpm.columns.map(sample_info.study),
        "cos_time": (numpy.cos(tpm.columns.map(sample_info.time) % 24) / 24 * 2 *numpy.pi),
        "sin_time": (numpy.sin(tpm.columns.map(sample_info.time) % 24) / 24 * 2 *numpy.pi),
    })
    fit = smf.mixedlm("tpm ~ seq_type + cos_time + sin_time", data=d, groups=d.study).fit()
