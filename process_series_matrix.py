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
    with gzip.open(series_matrix_path, "rt", newline='') as series_matrix_stream:
        series_matrix = series_matrix_stream.read()
        # Remove carriage returns
        # some files have these sprinkled throughout, in the middle of lines!
        series_matrix = series_matrix.replace('\r', '')

        # Read all the entries in the matrix
        # and extract the parts for the series and for the samples separately
        for line in series_matrix.splitlines():
            line = line.strip()
            if line.startswith("!Series_"):
                key, *values = line.removeprefix("!Series_").split("\t")
                series_data[key] = [strip_quotes(value) for value in values]
            if line.startswith("!Sample_"):
                key, *values = line.removeprefix("!Sample_").split("\t")
                # Some keys appear multiple times, so we will suffix later ones with 1, 2, etc
                key_count = len([k for k in sample_data.keys() if k.startswith(key)])
                if key_count == 0:
                    key_count = ''
                sample_data[key+f"{key_count}"] = [strip_quotes(value) for value in values]

    sample_data = pandas.DataFrame.from_dict(sample_data, 'columns', dtype="str")
    series_data = pandas.DataFrame.from_dict(series_data, 'index', dtype="str")
    print(f"Found {sample_data.shape[1]} values for {sample_data.shape[0]} samples")

    sample_data['SRX'] = ''
    for col in sample_data.columns:
        if col.startswith("relation"):
            srx = sample_data[col].map(extract_SRX)
            if any(srx != ''):
                sample_data['SRX'] = srx
    print(f"Identified SRX values for {sum(sample_data.SRX != '')} samples")

    return series_data, sample_data
