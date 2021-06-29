import pathlib
import pandas
from  process_series_matrix import process_series_matrix

targets = [
    {
        "GSE": "GSE165198",
        "sample_selector": lambda x: True,
    },
    {
        "GSE": "GSE70497",
        "sample_selector": lambda x: "WT" in x["source_name_ch1"],
    },
    {
        "GSE": "GSE40190",
        "sample_selector": lambda x: True,
    },
]
targets_by_gse = {record['GSE']:record for record in targets}

wildcard_constraints:
    GSE = "GSE(\d+)"

rule all:
    input:
        expand("data/{GSE}/sample_data.txt", GSE=targets_by_gse.keys()),
        #"data/GSE70497/salmon/GSM1787223/",
        "data/GSE40190/expression.txt"

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

checkpoint split_samples:
    input:
        "data/{GSE}/sample_data.txt"
    output:
        directory("data/{GSE}/samples")
    run:
        outdir = pathlib.Path(output[0])
        outdir.mkdir(exist_ok=True)

        sample_data = pandas.read_csv(input[0], sep="\t")
        print("Sample data", sample_data.source_name_ch1)
        selector = targets_by_gse[wildcards.GSE]['sample_selector']
        SRX = {row['geo_accession']:row['SRX'] for i, row in sample_data.iterrows() if selector(row)}
        for sample, srx in SRX.items():
            print(f"Processing {sample}, {srx}")
            sample_dir = (outdir / sample)
            sample_dir.mkdir(exist_ok=True)
            (sample_dir / "SRX.txt").open("w").write(srx)

rule download_select_sra_files:
    input:
        "data/{GSE}/samples/{sample}/SRX.txt"
    output:
        directory("data/{GSE}/SRA/{sample}/")
    message:
        "Fetching SRA files for {wildcards.GSE}"
    run:
        pathlib.Path(output[0]).mkdir(exist_ok=True)
        with pathlib.Path(input[0]).open() as srx_file:
            srx = srx_file.read().strip()
        if srx != '':
            shell(f"prefetch -O data/{wildcards.GSE}/SRA/{wildcards.sample}/ {srx}")

rule extract_fastq:
    input:
        "data/{GSE}/SRA/{sample}"
    output:
        directory("data/{GSE}/fastq/{sample}")
    message:
        "Extracting SRA to FastQ for {wildcards.GSE} {wildcards.sample}"
    shell:
        "fastq-dump --split-files -O data/{wildcards.GSE}/fastq/{wildcards.sample} {input}/*/*.sra"

rule run_salmon:
    input:
        "data/{GSE}/fastq/{sample}"
    output:
        directory("data/{GSE}/salmon/{sample}")
    resources:
        mem_mb=20000
        threads=6
    shell:
        "salmon quant -l IU -i /project/itmatlab/for_dimitra/pseudoalign_benchmark/dimitra/RevisionBMC/annotation/salmon.index/Mus_musculus.GRCm38.75 -1 {input}/*_1.fastq -2 {input}/*_2.fastq -o {output} -p 6"

def all_selected_samples(GSE):
    output = checkpoints.split_samples.get(GSE=GSE).output[0]
    samples_dir = pathlib.Path(output)
    return [sample_dir.name for sample_dir in samples_dir.glob("GSM*")]

def all_salmon_output(wildcards):
    samples = all_selected_samples(wildcards.GSE)
    return [f"data/{wildcards.GSE}/salmon/{sample}" for sample in samples]

rule aggregate_expression_values:
    input:
        all_salmon_output
    output:
        "data/{GSE}/expression.txt"
    run:
        pass
