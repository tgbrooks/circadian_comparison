import pathlib
import json
import re
import pandas
from  scripts.process_series_matrix import process_series_matrix

targets = {
    "Schwartz21": {
        "GSE": "GSE165198",
        "sample_selector": lambda x: True,
    },

    "Yang16A": {
        "GSE": "GSE70497",
        "sample_selector": lambda x: x.genotype == "WT",
    },

    "Lahens15": {
        "GSE": "GSE40190",
        "sample_selector": lambda x: True,
    },

    "Weger18": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x['gut microbiome status'] == "conventional raised"
    },
}

wildcard_constraints:
    sample = "GSM(\d+)"

rule all:
    input:
        expand("data/{study}/sample_data.txt", study=targets.keys()),
        "data/Lahens15/label_expression.tpm.txt",
        "data/Lahens15/label_expression.num_reads.txt",
        "data/Lahens15/salmon.meta_info.json",
        "data/Lahens15/jtk/JTKresult_expression.tpm.txt",
        "data/Weger18/label_expression.tpm.txt",
        "data/Weger18/label_expression.num_reads.txt",
        "data/Weger18/salmon.meta_info.json",
        "data/Weger18/jtk/JTKresult_expression.tpm.txt",

rule get_series_matrix:
    output:
        "data/{study}/series_matrix.txt.gz"
    message:
        "Fetching Series Matrix for {wildcards.study}"
    run:
        # Like GSE1nnn, the folder where the GSE series matrix resides
        GSE = targets[wildcards.study]['GSE']
        gse_folder = GSE[:-3]+"nnn"
        shell(f"wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/{gse_folder}/{GSE}/matrix/{GSE}_series_matrix.txt.gz -O {output}")

rule process_series_matrix:
    input:
        "data/{study}/series_matrix.txt.gz"
    output:
        "data/{study}/series_data.txt",
        "data/{study}/sample_data.txt"
    message:
        "Processing series matrix for {wildcards.study}"
    run:
        series_data, sample_data = process_series_matrix(input[0])
        series_data.to_csv(output[0], sep="\t", header=False)
        sample_data.to_csv(output[1], sep="\t", index=False)

checkpoint split_samples:
    input:
        "data/{study}/sample_data.txt"
    output:
        directory("data/{study}/samples")
    run:
        outdir = pathlib.Path(output[0])
        outdir.mkdir(exist_ok=True)

        sample_data = pandas.read_csv(input[0], sep="\t")
        print("Sample data", sample_data.source_name_ch1)
        selector = targets[wildcards.study]['sample_selector']
        SRX = {row['geo_accession']:row['SRX'] for i, row in sample_data.iterrows() if selector(row)}
        for sample, srx in SRX.items():
            print(f"Processing {sample}, {srx}")
            sample_dir = (outdir / sample)
            sample_dir.mkdir(exist_ok=True)
            (sample_dir / "SRX.txt").open("w").write(srx)

rule download_sra_files:
    input:
        "data/{study}/samples/{sample}/SRX.txt"
    output:
        directory("data/{study}/SRA/{sample}/")
    message:
        "Fetching SRA files for {wildcards.study}"
    run:
        pathlib.Path(output[0]).mkdir(exist_ok=True)
        with pathlib.Path(input[0]).open() as srx_file:
            srx = srx_file.read().strip()
        if srx != '':
            shell(f"prefetch -O data/{wildcards.study}/SRA/{wildcards.sample}/ {srx}")

rule extract_fastq:
    input:
        "data/{study}/SRA/{sample}"
    output:
        directory("data/{study}/fastq/{sample}")
    message:
        "Extracting SRA to FastQ for {wildcards.study} {wildcards.sample}"
    shell:
        "fastq-dump --split-files -O data/{wildcards.study}/fastq/{wildcards.sample} {input}/*/*.sra"

rule generate_salmon_index:
    output:
        directory("index/mouse_k{k}")
    resources:
        mem_mb=70000,
        threads=16
    shell:
        "salmon index -p 16 -i {output} \
             -t /project/itmatlab/for_dimitra/pseudoalign_benchmark/dimitra/RevisionBMC/annotation/salmon.index/Mus_musculus.GRCm38.75.TranscriptSeq.std.merged_with.dna.primary_assembly.fa \
             -d /project/itmatlab/for_dimitra/pseudoalign_benchmark/dimitra/RevisionBMC/annotation/salmon.index/decoy_names.txt \
             -k {wildcards.k}"

rule run_salmon:
    input:
        "data/{study}/fastq/{sample}",
        "index/mouse_k31"
    output:
        directory("data/{study}/salmon/{sample}")
    params:
        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",
        args = "-l A --softclip --softclipOverhangs --seqBias --gcBias --reduceGCMemory --biasSpeedSamp 10 --posBias -p 6"
    message:
        "Salmon: quantify {wildcards.sample} from {wildcards.study}"
    resources:
        mem_mb=25000,
        threads=6
    shell:
        "salmon quant -i {input[1]} -g {params.gtf_file} {params.args} -1 {input[0]}/*_1.fastq -2 {input[0]}/*_2.fastq -o {output}"

def all_selected_samples(study):
    ''' List all selected sample identifiers for a study '''
    output = checkpoints.split_samples.get(study=study).output[0]
    samples_dir = pathlib.Path(output)
    return [sample_dir.name for sample_dir in samples_dir.glob("GSM*")]

def all_salmon_output(wildcards):
    ''' List the salmon output files from all (selected) samples '''
    samples = all_selected_samples(wildcards.study)
    return [f"data/{wildcards.study}/salmon/{sample}" for sample in samples]

rule aggregate_expression_values:
    input:
        all_salmon_output
    output:
        "data/{study}/expression.tpm.txt",
        "data/{study}/expression.num_reads.txt",
        "data/{study}/salmon.meta_info.json"
    message:
        "Aggregate Salmon quantifications for {wildcards.study}"
    run:
        tpm_dict = {}
        num_reads_dict = {}
        meta_info = {}
        for sample,sampledir in zip(all_selected_samples(wildcards.study),input):
            samplequantfile = sampledir + "/quant.genes.sf"
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

rule label_data:
    input:
        "data/{study}/expression.tpm.txt",
        "data/{study}/expression.num_reads.txt",
        "data/{study}/sample_data.txt"
    output:
        "data/{study}/label_expression.tpm.txt",
        "data/{study}/label_expression.num_reads.txt"
    run:
        sample = pandas.read_csv(input[2], sep="\t", index_col="geo_accession")
        tpm = pandas.read_csv(input[0], sep="\t")
        tpm.columns = sample.reindex(tpm.columns, fill_value="Name").title
        tpm.to_csv(output[0], sep="\t", index=False)
        num_reads = pandas.read_csv(input[1], sep="\t")
        num_reads.columns = sample.reindex(num_reads.columns, fill_value="Name").title
        num_reads.to_csv(output[1], sep ="\t", index=False)

def sample_timepoints(study):
    sample_data = pandas.read_csv(f"data/{study}/sample_data.txt", sep="\t", index_col="geo_accession")
    expression_table = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
    times = list(sample_data.loc[expression_table.columns].time)
    times = [int(re.match("[ZC]T(\d+)", time).groups()[0]) for time in times]
    return times

rule run_JTK:
    input:
        "data/{study}/expression.tpm.txt",
        "data/{study}/sample_data.txt"
    output:
        "data/{study}/jtk/JTKresult_expression.tpm.txt"
    params:
        timepoints = lambda wildcards: sample_timepoints(wildcards.study),
        out_dir = "data/{study}/jtk/"
    script:
        "scripts/run_jtk.R"
