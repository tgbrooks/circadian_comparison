'''
Perform outlier detection.

Outliers are found by PCA within each study.
Samples that lie > 3 standard deviations from the mean according
to the first principal component are dropped.
'''

import pandas
import numpy
import scipy.sparse.linalg

ZSCORE_OUTLIER_THRESHOLD = 3

data = pandas.read_csv(snakemake.input.tpm, sep="\t", index_col=[0,1])
sample_info = pandas.read_csv(snakemake.input.sample_info, sep="\t", index_col=[0])
print(data.shape)
print(sample_info.shape)
assert (data.columns == sample_info.index).all()

data = numpy.log(data + 0.01) # Log transform the data

def standardize(x):
    ''' z-score data '''
    std =  x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
    std[std.isna().any(axis=1)] = 0 # Constant rows should give 0s not NaNs
    std[x.std(axis=1) < 1e-10] = 0
    return std

outliers = set()
for study, study_data in data.groupby(sample_info.study, axis=1):

    is_outlier = pandas.Series(False, index=study_data.columns)
    while True:
        drop = is_outlier[study_data.columns]
        std_data = standardize(study_data.loc[:,~drop])
        u, s, dt = scipy.sparse.linalg.svds(std_data, k=1)

        scores = (std_data.T @ u)[0]
        scores_std = numpy.std(scores)
        is_new_outlier = scores.abs()  - scores.mean() > ZSCORE_OUTLIER_THRESHOLD * scores_std
        if sum(is_new_outlier) == 0:
            break
        is_outlier[is_new_outlier.index[is_new_outlier]] = True
    outliers = outliers.union(is_outlier.index[is_outlier])
    print(f"Identified {sum(is_outlier)} in {study}")
    print('\n'.join(is_outlier.index[is_outlier]))

print(f"Identified {len(outliers)} in total")
pandas.Series(list(outliers)).to_csv(snakemake.output.outliers, sep="\t", index=False, header=False)
