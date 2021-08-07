import pathlib
import json
import re
import pandas
from  scripts.process_series_matrix import process_series_matrix

from studies import targets, studies, sample_timepoints, tissues

wildcard_constraints:
    sample = "GSM(\d+)",

rule all:
    input:
        # All study-level files:
        expand("data/{study}/sample_data.txt", study=targets.keys()),
        expand("data/{study}/label_expression.tpm.txt", study=studies),
        expand("data/{study}/jtk/JTKresult_expression.tpm.txt", study=studies),
        # All tissue-level files:
        expand("results/{tissue}/{file}",
            tissue = tissues,
            file = [
                "qc/percent_mapping.png",
                "clock_genes/plot_Arntl.png",
                "jtk/breakdowns.png",
                "num_common_genes.txt",
                "PCA/all_samples_study.png",
                "robustness/expression_level_robust.png",
                "tpm_all_samples.txt",
                "amplitude_scatter_grid.png",
            ]
        )

rule get_series_matrix:
    output:
        "data/{study}/series_matrix.txt.gz"
    message:
        "Fetching Series Matrix for {wildcards.study}"
    run:
        # Like GSE1nnn, the folder where the GSE series matrix resides
        GSE = targets[wildcards.study]['GSE']
        gse_folder = GSE[:-3]+"nnn"
        postfix = targets[wildcards.study].get('GSE_postfix', '')
        shell(f"wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/{gse_folder}/{GSE}/matrix/{GSE}{postfix}_series_matrix.txt.gz -O {output}")

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
        if len(SRX) == 0:
            raise Exception(f"No samples selected for {wildcards.study}")
        for sample, srx in SRX.items():
            print(f"Processing {sample}, {srx}")
            sample_dir = (outdir / sample)
            sample_dir.mkdir(exist_ok=True)
            (sample_dir / "SRX.txt").open("w").write(srx)

rule download_sra_files:
    input:
        "data/{study}/samples/{sample}/SRX.txt"
    output:
        temp(directory("data/{study}/SRA/{sample}/"))
    message:
        "Fetching SRA files for {wildcards.study}:{wildcards.sample}"
    resources:
        ncbi_download=1,
    run:
        pathlib.Path(output[0]).mkdir(exist_ok=True)
        srrfile = f"data/{wildcards.study}/samples/{wildcards.sample}/SRR.txt"
        with pathlib.Path(input[0]).open() as srx_file:
            srx = srx_file.read().strip()
        if srx != '':
            # First get the SRR numbers from the SRX
            shell(f'efetch -db sra -format runinfo -id {srx} | grep -oh "^SRR[0-9]\\+" > {srrfile}')
            shell(f"prefetch -O data/{wildcards.study}/SRA/{wildcards.sample}/ `cat {srrfile}`")

rule extract_fastq:
    input:
        "data/{study}/SRA/{sample}"
    output:
        temp(directory("data/{study}/fastq/{sample}"))
    priority: 10
    message:
        "Extracting SRA to FastQ for {wildcards.study} {wildcards.sample}"
    shell:
        # Note: Refer to https://edwards.sdsu.edu/research/fastq-dump/ for information about using fastq-dump properly
        "fastq-dump --readids --skip-technical --split-files --clip -F -I -O data/{wildcards.study}/fastq/{wildcards.sample} {input}/*/*.sra"

rule generate_salmon_index:
    output:
        directory("index/mouse_k{k}")
    threads: 16
    resources:
        mem_mb=70000,
    shell:
        "salmon index -p 16 -i {output} \
             -t /project/itmatlab/index/SALMON-1.4.0/Mus_musculu.GRCm38.75/Mus_musculus.GRCm38.75.TranscriptSeq.std.merged_with.dna.primary_assembly.fa \
             -d /project/itmatlab/index/SALMON-1.4.0/Mus_musculu.GRCm38.75/index/decoy_names.txt \
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
    run:
        fastqdir=pathlib.Path(input[0])
        if len(list(fastqdir.glob("*_2.fastq")))>0:
            shell(f"salmon quant -i {input[1]} -g {params.gtf_file} {params.args} -1 {input[0]}/*_1.fastq -2 {input[0]}/*_2.fastq -o {output[1]}")
        else:
            shell(f"salmon quant -i {input[1]} -g {params.gtf_file} {params.args} -r {input[0]}/*_1.fastq -o {output[1]}")

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
        expression_tpm = expand("data/{study}/label_expression.tpm.txt", study=studies),
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
    params:
        studies = studies,
    output:
        percent_mapping = "results/{tissue}/qc/percent_mapping.png",
        num_processed = "results/{tissue}/qc/num_processed.png",
        num_mapped = "results/{tissue}/qc/num_mapped.png",
        expression_distributions = "results/{tissue}/qc/expression_distributions.png",
        ENSMUSG00000086503 = "results/{tissue}/qc/Xist_expression.png",
        ENSMUSG00000029368 = "results/{tissue}/qc/Alb_expression.png",
        ENSMUSG00000064337  = "results/{tissue}/qc/mtRnr1_expression.png",
        ENSMUSG00000064339 = "results/{tissue}/qc/mtRnr2_expression.png",
        ENSMUSG00000064336 = "results/{tissue}/qc/mtTf_expression.png",
    resources:
        mem_mb = 4000,
    script:
        "scripts/qc_plots.py"

rule plot_genes:
    input:
        expression_tpm = expand("data/{study}/label_expression.tpm.txt", study=studies),
    params:
        studies = studies,
        genes_ID = ["ENSMUSG00000055116", "ENSMUSG00000020038", "ENSMUSG00000068742", "ENSMUSG00000020893", "ENSMUSG00000055866", "ENSMUSG00000028957", "ENSMUSG00000020889", "ENSMUSG00000021775", "ENSMUSG00000059824", "ENSMUSG00000029238"],
        genes_symbol = ["Arntl", "Cry1", "Cry2", "Per1", "Per2", "Per3", "Nr1d1", "Nr1d2", "Dbp", "Clock"]
    output:
        ENSMUSG00000055116 = "results/{tissue}/clock_genes/plot_Arntl.png",
        ENSMUSG00000020038 = "results/{tissue}/clock_genes/plot_Cry1.png",
        ENSMUSG00000068742 = "results/{tissue}/clock_genes/plot_Cry2.png",
        ENSMUSG00000020893 = "results/{tissue}/clock_genes/plot_Per1.png",
        ENSMUSG00000055866 = "results/{tissue}/clock_genes/plot_Per2.png",
        ENSMUSG00000028957 = "results/{tissue}/clock_genes/plot_Per3.png",
        ENSMUSG00000020889 = "results/{tissue}/clock_genes/plot_Nr1d1.png",
        ENSMUSG00000021775 = "results/{tissue}/clock_genes/plot_Nr1d2.png",
        ENSMUSG00000059824 = "results/{tissue}/clock_genes/plot_Dbp.png",
        ENSMUSG00000029238 = "results/{tissue}/clock_genes/plot_Clock.png",
        ENSMUSG00000055116mod24 = "results/{tissue}/genes24/plot_Arntl.png",
        ENSMUSG00000020038mod24 = "results/{tissue}/genes24/plot_Cry1.png",
        ENSMUSG00000068742mod24 = "results/{tissue}/genes24/plot_Cry2.png",
        ENSMUSG00000020893mod24 = "results/{tissue}/genes24/plot_Per1.png",
        ENSMUSG00000055866mod24 = "results/{tissue}/genes24/plot_Per2.png",
        ENSMUSG00000028957mod24 = "results/{tissue}/genes24/plot_Per3.png",
        ENSMUSG00000020889mod24 = "results/{tissue}/genes24/plot_Nr1d1.png",
        ENSMUSG00000021775mod24 = "results/{tissue}/genes24/plot_Nr1d2.png",
        ENSMUSG00000059824mod24 = "results/{tissue}/genes24/plot_Dbp.png",
        ENSMUSG00000029238mod24 = "results/{tissue}/genes24/plot_Clock.png",
    resources:
        mem_mb = 4000,
    script:
        "scripts/gene_plots.py"

rule plot_jtk:
    input:
        jtk = expand("data/{study}/jtk/JTKresult_expression.tpm.txt", study=studies),
    params:
        studies = studies,
    output:
        breakdowns = "results/{tissue}/jtk/breakdowns.png",
        periods = "results/{tissue}/jtk/periods.png",
        amplitudes = "results/{tissue}/jtk/amplitudes.png",
        phases = "results/{tissue}/jtk/phases.png",
    script:
        "scripts/plot_jtk.py"

rule plot_overlapped_genes:
    input:
        jtk = expand("data/{study}/jtk/JTKresult_expression.tpm.txt", study=studies),
        tpm = expand("data/{study}/expression.tpm.txt", study=studies),
    params:
        studies = studies,
    output:
        num_common_genes = "results/{tissue}/num_common_genes.txt",
        common_genes_pvalue = "results/{tissue}/common_genes_pvalue.txt",
        robustness_score = "results/{tissue}/robustness_score.txt",
        heatmap = "results/{tissue}/common_genes_heatmap.png",
    resources:
        mem_mb = 4000
    script:
        "scripts/plot_overlapped_genes.py"

rule plot_PCA:
    input:
        tpm = expand("data/{study}/expression.tpm.txt", study=studies),
    params:
        studies = studies,
    resources:
        mem_mb = 6000,
    output:
        all_samples_study = "results/{tissue}/PCA/all_samples_study.png",
        all_samples_time = "results/{tissue}/PCA/all_samples_time.png",
    script:
        "scripts/plot_PCA.py"

rule robustness:
    input:
        robustness = "results/robustness_score.txt",
        tpm = expand("data/{study}/expression.tpm.txt",study=studies),
        jtk = expand("data/{study}/jtk/JTKresult_expression.tpm.txt", study=studies),
    params:
        studies = studies,
    output:
        expression_level = "results/{tissue}/robustness/expression_level_robust.png",
        amplitude = "results/{tissue}/robustness/amplitude_robust.png",
        period = "results/{tissue}/robustness/period_robust.png",
        phase = "results/{tissue}/robustness/phase_robust.png",
    script:
        "scripts/plot_robustness.py"

rule all_samples:
    input:
        tpm = expand("data/{study}/expression.tpm.txt",study=studies),
        num_reads = expand("data/{study}/expression.num_reads.txt",study=studies),
    params:
        studies = studies,
    output:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        num_reads = "results/{tissue}/num_reads_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
    run:
        all_tpm = pandas.DataFrame()
        all_num_reads = pandas.DataFrame()
        all_samples = []
        for study, tpmfile in zip(studies, input.tpm):
            if study == "Weger18":
                continue
            tpm = pandas.read_csv(tpmfile, sep="\t", index_col="Name")
            time = sample_timepoints(study)
            all_samples.append(pandas.DataFrame({
            "study": [study for i in range(len(tpm.columns))],
            "time": time
            }, index=tpm.columns))
            all_tpm = pandas.concat([all_tpm, tpm], axis=1)
        all_tpm.insert(0, 'Symbol', pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol'])
        all_tpm.to_csv(output[0], sep ="\t", index="Name")
        all_samples_concat = pandas.concat(all_samples, axis=0)
        all_samples_df = pandas.DataFrame(all_samples_concat)
        all_samples_df.to_csv(output[2], sep="\t")

        for study, numreadsfile in zip(studies, input.num_reads):
            if study == "Weger18":
                continue
            num_reads = pandas.read_csv(numreadsfile, sep="\t", index_col="Name")
            all_num_reads = pandas.concat([all_num_reads, num_reads], axis=1)
        all_num_reads.insert(0, 'Symbol', pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol'])
        all_num_reads.to_csv(output[1], sep ="\t", index="Name")

rule scatter_grid:
    input:
        jtk = expand("data/{study}/jtk/JTKresult_expression.tpm.txt", study=studies),
    params:
        studies = studies,
    output:
        amplitude = "results/{tissue}/amplitude_scatter_grid.png",
        period = "results/{tissue}/period_scatter_grid.png",
        phase = "results/{tissue}/phase_scatter_grid.png",
    resources: 
        mem_mb = 12000
    script:
        "scripts/plot_scatter_grid.py"
