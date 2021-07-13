import pathlib
import json
import re
import pandas
from  scripts.process_series_matrix import process_series_matrix

from studies import targets, studies

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
        "data/Weger18_Liver_M/label_expression.tpm.txt",
        "data/Weger18_Liver_M/label_expression.num_reads.txt",
        "data/Weger18_Liver_M/salmon.meta_info.json",
        "data/Weger18_Liver_M/jtk/JTKresult_expression.tpm.txt",
        "data/Weger18_Liver_F/label_expression.tpm.txt",
        "data/Weger18_Liver_F/label_expression.num_reads.txt",
        "data/Weger18_Liver_F/salmon.meta_info.json",
        "data/Weger18_Liver_F/jtk/JTKresult_expression.tpm.txt",
        "data/Zhang14_RNAseq_Liver_M/label_expression.tpm.txt",
        "data/Zhang14_RNAseq_Liver_M/label_expression.num_reads.txt",
        "data/Zhang14_RNAseq_Liver_M/salmon.meta_info.json",
        "data/Zhang14_RNAseq_Liver_M/jtk/JTKresult_expression.tpm.txt",
        "results/qc.percent_mapping.png",

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
        # Note: Refer to https://edwards.sdsu.edu/research/fastq-dump/ for information about using fastq-dump properly
        "fastq-dump --readids --skip-technical --split-files --clip -O data/{wildcards.study}/fastq/{wildcards.sample} {input}/*/*.sra"

rule generate_salmon_index:
    output:
        directory("index/mouse_k{k}")
    threads: 16
    resources:
        mem_mb=70000,
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
        "data/{study}/salmon/{sample}/quant.genes.sf",
        directory("data/{study}/salmon/{sample}")
    params:
        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",
        args = "-l A --softclip --softclipOverhangs --seqBias --gcBias --reduceGCMemory --biasSpeedSamp 10 --posBias -p 6"
    message:
        "Salmon: quantify {wildcards.sample} from {wildcards.study}"
    threads: 6
    resources:
        mem_mb=25000,
    shell:
        "salmon quant -i {input[1]} -g {params.gtf_file} {params.args} -1 {input[0]}/*_1.fastq -2 {input[0]}/*_2.fastq -o {output[1]}"

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

rule extract_GeneSymbol:
    input:
        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf"
    output:
        "gene_name.txt"
    run:
        shell("""cat {input} | awk 'BEGIN{{FS="\\t"}}{{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]}}' | sed 's/gene_id "//g' | sed 's/gene_name "//g' | sed 's/"//g' | sed 's/ //g' > {output}""")
        id_name = pandas.read_csv("gene_name.txt", sep="\t", header=None)
        id_name.to_csv("gene_name.txt", sep="\t",header=["ID","GeneSymbol"])

rule label_data:
    input:
        "data/{study}/expression.tpm.txt",
        "data/{study}/expression.num_reads.txt",
        "data/{study}/sample_data.txt",
        "gene_name.txt"
    output:
        "data/{study}/label_expression.tpm.txt",
        "data/{study}/label_expression.num_reads.txt"
    run:
        sample = pandas.read_csv(input[2], sep="\t", index_col="geo_accession")
        tpm = pandas.read_csv(input[0], sep="\t")
        tpm.columns = sample.reindex(tpm.columns, fill_value="ID").title
        gene_name_from_id = pandas.read_csv(input[3], sep="\t", index_col="ID")['GeneSymbol']
        tpm.insert(1, 'GeneSymbol', tpm['ID'].map(gene_name_from_id))
        tpm.to_csv(output[0], sep="\t", index=False)
        num_reads = pandas.read_csv(input[1], sep="\t")
        num_reads.columns = sample.reindex(num_reads.columns, fill_value="ID").title
        gene_name_from_id = pandas.read_csv(input[3], sep="\t", index_col="ID")['GeneSymbol']
        num_reads.insert(1, 'GeneSymbol', num_reads['ID'].map(gene_name_from_id))
        num_reads.to_csv(output[1], sep ="\t", index=False)

def sample_timepoints(study):
    sample_data = pandas.read_csv(f"data/{study}/sample_data.txt", sep="\t", index_col="geo_accession")
    expression_table = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
    times = targets[study]["time"](sample_data, expression_table)
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

rule plot_qc:
    input:
        salmon_metainfo = expand("data/{study}/salmon.meta_info.json", study=studies),
    params:
        studies = studies,
    output:
        percent_mapping = "results/qc.percent_mapping.png",
        num_processed = "results/qc.num_processed.png",
        num_mapped = "results/qc.num_mapped.png",
    script:
        "scripts/qc_plots.py"

rule plot_arntl:
    input:
        expression_tpm = expand("data/{study}/label_expression.tpm.txt", study=studies),
    params:
        studies = studies,
    output:
        arntl_expression_tpm = "results/plot_arntl.png",
    script:
        "scripts/arntl_plots.py"
