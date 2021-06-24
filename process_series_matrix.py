import gzip
import re

import pandas

def strip_quotes(string):
    return string.removeprefix('"').removesuffix('"')

def extract_SRX(string):
    match = re.search(r"(SRX\d+)", string)
    if match:
        return match.groups()[0]
    return ''

def process_series_matrix(series_matrix_path):
    series_data = {}
    sample_data = {}
    with gzip.open(series_matrix_path, "rt") as series_matrix:
        for line in series_matrix.readlines():
            line = line.strip()
            if line.startswith("!Series_"):
                key, *values = line.removeprefix("!Series_").split("\t")
                series_data[key] = [strip_quotes(value) for value in values]
            if line.startswith("!Sample_"):
                key, *values = line.removeprefix("!Sample_").split("\t")
                sample_data[key] = [strip_quotes(value) for value in values]

    sample_data = pandas.DataFrame.from_dict(sample_data, 'columns', dtype="str")
    series_data = pandas.DataFrame.from_dict(series_data, 'index', dtype="str")

    if "relation" in sample_data:
        sample_data['SRX'] = sample_data.relation.map(extract_SRX)
    else:
        sample_data['SRX'] = ''

    return series_data, sample_data
