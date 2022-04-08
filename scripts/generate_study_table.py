# Generate a table of the study info from the 'targets' object
import re
import json
import pandas
from studies import targets
studies = snakemake.params.studies
sample_info = pandas.read_csv(snakemake.input.sample_info, sep="\t", index_col=0)
outlier_samples = [x.strip() for x in open(snakemake.input.outliers).readlines()]
j = open(snakemake.input.study_classification).read()
study_seq_class = json.loads(j)

series_data = { study: open(f"data/{study}/series_data.txt").read() for study in studies }
def get_pubmed_ids(series_data):
    m = re.findall("pubmed_id\t([0-9]+)", series_data)
    return m
pubmed_ids = {study: get_pubmed_ids(series_data[study]) for study in studies}
sample_data = {study: open(f"data/{study}/sample_data.txt").read() for study in studies}
def sequence_type(study):
    study_samples = sample_info[sample_info.study == study].index
    study_samples = [line for line in sample_data[study].splitlines()
                        if line.split('\t')[1] in study_samples]
    value = []
    if any(re.search("TruSeq|TrueSeq|Tru-Seq", line, re.I) for line in study_samples):
        value.append('TruSeq')
    if any(re.search("RiboZero", line, re.I) for line in study_samples):
        value.append('RiboZero')
    if any(re.search("Total RNA", line, re.I) for line in study_samples):
        value.append('Total')
    if any(re.search("Stranded|stranded", line, re.I) for line in study_samples):
        value.append('Stranded')
    return ' '.join(value)
inferred_sequencing = {study: sequence_type(study) for study in studies}

def full_seq(study):
    seq = targets[study]['seq']
    seq_class = study_seq_class[study]
    stranded = 'US' if 'unstranded' in seq_class else 'SS'
    ends = 'SE' if 'unpaired' in seq_class else 'PE'
    full_seq = f"{stranded} {ends} {seq}"
    return full_seq

samples = {study: sample_info[sample_info.study == study] for study in studies}
def format_range(low, high):
    if not low == low: # NaN
        return 'unknown'
    if low == high:
        return str(low)
    else:
        return f"{low} - {high}"

sample_table = pandas.DataFrame({
    "Name": {study: targets[study]['short_name'] for study in studies},
    "GSE": {study: targets[study]['GSE'] for study in studies},
    "PUBMED": {study: ','.join(pubmed_ids[study]) for study in studies},
    "Sex": {study: targets[study]['sex'] for study in studies},
    "Light": {study: targets[study]['light'] for study in studies},
    "Age (weeks)": {
        study: format_range(targets[study]['age_low'], targets[study]['age_high'])
            for study in studies
    },
    "Sample Count": {
        study: len(samples[study])
            for study in studies
    },
    "Outliers": {
        study: sum(samples[study].index.isin(outlier_samples))
            for study in studies
    },
    "Timepoints Per Cycle": {
        study: (samples[study].time % 24).nunique()
            for study in studies
    },
    "Replicates": {
        study: format_range((samples[study].time).value_counts().min(), (samples[study].time).value_counts().max())
            for study in studies
    },
    "Cycles": {
        study: (samples[study].time).nunique() / (samples[study].time % 24).nunique()
            for study in studies
    },
    "Sequencing Type": {study: full_seq(study) for study in studies},
    "Inferred Sequencing Type": inferred_sequencing,
    "Note": {study: targets[study].get('note', '') for study in studies},
})

sample_table.sort_values(by="Name").to_csv(snakemake.output.table, sep="\t", index=False)
