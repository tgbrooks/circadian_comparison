import pathlib
import json
import pandas
from  process_series_matrix import process_series_matrix

targets = [
    {
        "GSE": "GSE165198",
        "sample_selector": lambda x: True,
    },
    {
        "GSE": "GSE70497",
        "sample_selector": lambda x: x.genotype == "WT",
    },
    {
        "GSE": "GSE40190",
        "sample_selector": lambda x: True,
    },
]
targets_by_gse = {record['GSE']:record for record in targets}

wildcard_constraints:
    GSE = "GSE(\d+)",
    sample = "GSM(\d+)"

rule all:
    input:
        expand("data/{GSE}/sample_data.txt", GSE=targets_by_gse.keys()),
        #"data/GSE70497/salmon/GSM1787223/",
        "data/GSE40190/expression.tpm.txt",
        "data/GSE40190/expression.num_reads.txt",
        "data/GSE40190/salmon.meta_info.json"

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

rule download_sra_files:
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
    message:
        "Salmon: quantify {wildcards.sample} from {wildcards.GSE}"
    resources:
        mem_mb=20000,
        threads=6
    shell:
        "salmon quant -l A -i /project/itmatlab/for_dimitra/pseudoalign_benchmark/dimitra/RevisionBMC/annotation/salmon.index/Mus_musculus.GRCm38.75 -g  /project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf -1 {input}/*_1.fastq -2 {input}/*_2.fastq -o {output} -p 6"

def all_selected_samples(GSE):
    ''' List all selected sample identifiers for a study '''
    output = checkpoints.split_samples.get(GSE=GSE).output[0]
    samples_dir = pathlib.Path(output)
    return [sample_dir.name for sample_dir in samples_dir.glob("GSM*")]

def all_salmon_output(wildcards):
    ''' List the salmon output files from all (selected) samples '''
    samples = all_selected_samples(wildcards.GSE)
    return [f"data/{wildcards.GSE}/salmon/{sample}" for sample in samples]

rule aggregate_expression_values:
    input:
        all_salmon_output
    output:
        "data/{GSE}/expression.tpm.txt",
        "data/{GSE}/expression.num_reads.txt",
        "data/{GSE}/salmon.meta_info.json"
    message:
        "Aggregate Salmon quantifications for {wildcards.GSE}"
    run:
        tpm_dict = {}
        num_reads_dict = {}
        meta_info = {}
        for sample,sampledir in zip(all_selected_samples(wildcards.GSE),input): 
            samplequantfile = sampledir + "/quant.sf"
            quant = pandas.read_csv(samplequantfile,sep="\t", index_col=0)
            tpm_dict[sample]=quant.TPM
            num_reads_dict[sample]=quant.NumReads
            with open(sampledir + "/aux_info/meta_info.json") as metainfofile:
                meta_info[sample] = json.load(metainfofile)
        tpm = pandas.DataFrame.from_dict(tpm_dict,orient="columns")
        tpm.to_csv(output[0],sep = "\t")
        num_reads = pandas.DataFrame.from_dict(num_reads_dict,orient="columns")
        num_reads.to_csv(output[1],sep = "\t")
        with open(output[2], "wt") as meta_info_out:
            json.dump(meta_info, meta_info_out, indent=4)
