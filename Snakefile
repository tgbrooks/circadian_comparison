import pathlib
from  process_series_matrix import process_series_matrix

rule all:
    input:
        #"data/GSE165198/GSE165198_series_matrix.txt.gz",
        "data/GSE165198/sample_data.txt",
        "data/GSE70497/sample_data.txt"

rule get_series_matrix:
    output:
        "data/{GSE}/{GSE}_series_matrix.txt.gz"
    message:
        "Fetching Series Matrix for {wildcards.GSE}"
    run:
        # Like GSE1nnn, the folder where the GSE series matrix resides
        gse_folder = wildcards.GSE[:-3]+"nnn"
        shell("wget -P data/{wildcards.GSE}/ ftp://ftp.ncbi.nlm.nih.gov/geo/series/{gse_folder}/{wildcards.GSE}/matrix/{wildcards.GSE}_series_matrix.txt.gz")

rule process_series_matrix:
    input:
        "data/{GSE}/{GSE}_series_matrix.txt.gz"
    output:
        "data/{GSE}/series_data.txt",
        "data/{GSE}/sample_data.txt"
    message:
        "Processing series matrix for {wildcards.GSE}"
    run:
        series_data, sample_data = process_series_matrix(input[0])
        series_data.to_csv(output[0], sep="\t", header=False)
        sample_data.to_csv(output[1], sep="\t", index=False)
