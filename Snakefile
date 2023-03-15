import pathlib
import json
import re
import pandas
import collections
from  scripts.process_series_matrix import process_series_matrix

import studies as studies_config
from studies import targets, studies, sample_timepoints, select_tissue, studies_by_tissue
def sample_timepoints(study, **kwargs):
    # Force checkpoint checking...
    output = checkpoints.split_samples.get(study=study).output[0]
    output = checkpoints.outlier_detection.get(study=study).output[0]
    return studies_config.sample_timepoints(study, **kwargs)

tissues = ["Liver"]

SPLINE_FIT_N_BATCHES = 400
# Require at least this number of mean reads per gene to include in q-value computations
MEAN_READCOUNT_THRESHOLD = 2
WORKING_DIR = "/project/itmatlab/for_tom/circadian_controls/"
BOOTEJTK_SIF = "~/.apptainer/images/ejtk_bootejtk.sif"

wildcard_constraints:
    sample = "GSM(\d+)",
    period = ".*", # Allow empty periods, empty means default
    permutation = r"(_perm/[0-9]+|)", # Permutation '' means no permutation, else a number
    study = "[a-zA-Z0-9_-]+",

rule all:
    input:
        # All study-level files:
        expand("data/{study}/sample_data.txt", study=targets.keys()),
        expand("data/{study}/label_expression.tpm.txt", study=studies),
        expand("data/{study}/jtk.results.txt", study=studies),
        expand("data/{study}/bootejtk/results.txt", study=studies),
        # All tissue-level files:
        expand("results/{tissue}/{file}",
            tissue = tissues,
            file = [
                "qc/percent_mapping.png",
                "clock_genes/plot_Arntl.png",
                "jtk/breakdowns.png",
                "jtk24/breakdowns.png",
                "jtk12/breakdowns.png",
                "jtk8/breakdowns.png",
                "jtk/num_common_genes.txt",
                "jtk24/num_common_genes.txt",
                "jtk12/num_common_genes.txt",
                "jtk8/num_common_genes.txt",
                "common_nonrhythmic_genes.txt",
                "PCA",
                "jtk24/robustness/expression_level_robust.png",
                "jtk12/robustness/expression_level_robust.png",
                "jtk8/robustness/expression_level_robust.png",
                "tpm_all_samples.txt",
                "jtk.results.txt",
                "jtk24.results.txt",
                "jtk12.results.txt",
                "jtk8.results.txt",
                #"amplitude_scatter_grid.png",
                "consensus_pca/",
                "outlier_samples.txt",
                "assess_jtk/period_statistics.txt",
                "study_table.txt",
            ]
        ),
        "results/Liver/compareRhythms/plots/",
        "results/Liver/compareRhythms/summary.txt",
        # NOTE: big computation, ~500 hours of CPU time
        "results/Liver/spline_fit/summary.txt",
        #"results/Liver/spline_fit_perm/1/batches/1.summary.txt",
        "results/Liver/spline_fit/tsne.png",
        "results/Liver/spline_fit/stats/",
        "results/Liver/spline_fit_perm/1/summary.txt",
        "results/Liver/spline_fit_perm/1/stats/",
        "results/Liver/spline_fit/phase_variability/phase_std_distribution.png",
        "results/Liver/stable_genes/stable_gene_list.txt",
        "results/Liver/supplemental/",

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
        "data/{study}/samples/split_flag"
    run:
        outflag = pathlib.Path(output[0])
        outdir = outflag.parent
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
        outflag.touch()

rule download_sra_files:
    input:
        "data/{study}/samples/{sample}/SRX.txt",
        "data/{study}/samples/split_flag",
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

rule trim_adapators:
    input:
        "data/{study}/fastq/{sample}/",
    output:
        temp(directory("data/{study}/fastq_trimmed/{sample}")),
    params:
        args = lambda wildcards: targets[wildcards.study]['trim_adaptor'],
    message:
        "Trimming FASTQ for {wildcards.study} {wildcards.sample}"
    script:
        "scripts/trim_adaptors.py"

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
        lambda wildcards: "data/{study}/fastq_trimmed/{sample}" if targets[wildcards.study].get("trim_adaptor", False) else "data/{study}/fastq/{sample}",
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
    samples_dir = pathlib.Path(output).parent
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

rule classify_studies:
    input:
        meta_info = expand("data/{study}/salmon.meta_info.json", study=studies)
    output:
        "results/study_classification.json"
    script:
        "scripts/classify_studies.py"

rule prep_jtk:
    input:
        tpm = "data/{study}/expression.tpm.txt",
        sample_data = "data/{study}/sample_data.txt",
        outliers = "data/{study}/outlier_samples.txt",
    output:
        tpm = temp("data/{study}/expression.tpm.for_JTK.txt")
    run:
        # Select just the non-outlier samples and put in a file so JTK can read it
        outlier_samples = [x.strip() for x in open(input.outliers).readlines()]
        tpm = pandas.read_csv(input.tpm, sep="\t", index_col=0)
        samples = tpm.columns
        tpm_selected = tpm[[sample for sample in samples if sample not in outlier_samples]]
        tpm_selected.to_csv(output.tpm, sep="\t")

rule run_jtk:
    input:
        "data/{study}/expression.tpm.for_JTK.txt",
        "data/{study}/sample_data.txt",
        "data/{study}/expression.tpm.txt",
        lambda wildcards: checkpoints.split_samples.get(study=wildcards.study).output[0], # Trigger checkpoint
        outliers = "data/{study}/outlier_samples.txt",
    output:
        "data/{study}/jtk{period}/JTKresult_expression.tpm.for_JTK.txt"
    params:
        timepoints = lambda wildcards: sample_timepoints(wildcards.study, drop_outliers=True),
        out_dir = "data/{study}/jtk{period}/",
        period = lambda wildcards: 'default' if wildcards.period == '' else wildcards.period
    script:
        "scripts/run_jtk.R"

rule process_jtk:
    input:
        "data/{study}/jtk{period}/JTKresult_expression.tpm.for_JTK.txt",
        "data/{study}/expression.num_reads.txt",
    output:
        "data/{study}/jtk{period}.results.txt",
    run:
        # This generates a JTK output with q-values computed after removing
        # low-expressed genes.
        jtk = pandas.read_csv(input[0], sep="\t", index_col=0)
        num_reads = pandas.read_csv(input[1], sep="\t", index_col=0).loc[jtk.index]
        selected = num_reads.mean(axis=1) >= MEAN_READCOUNT_THRESHOLD
        print(num_reads)
        print(jtk)
        ps = jtk.loc[selected, 'ADJ.P']
        import statsmodels.api as sm
        print(selected)
        print(ps)
        if len(ps) > 0:
            _, qs, _, _ = sm.stats.multipletests(ps, method="fdr_bh")
        else:
            ps = []
        jtk['dropped'] = ~selected# Mark any genes as 'dropped' if they did not pass the cutoff
        jtk['qvalue'] = 1 # Dropped genes get q=1
        jtk.loc[selected, 'qvalue'] = qs
        jtk.to_csv(output[0], sep="\t")

rule prep_bootejtk:
    input:
        tpm = "data/{study}/expression.tpm.txt",
        sample_data = "data/{study}/sample_data.txt",
        outliers = lambda wildcards: f"results/{targets[wildcards.study]['tissue']}/outlier_samples.txt",
        split_samples = lambda wildcards: checkpoints.split_samples.get(study=wildcards.study).output[0], # Trigger checkpoint
    output:
        tpm = temp("data/{study}/bootejtk/expression.tpm.for_BooteJTK.txt",)
    params:
        timepoints = lambda wildcards: sample_timepoints(wildcards.study, drop_outliers=True),
    run:
        # Select just the non-outlier samples and put in a file so BooteJTK can read it
        outlier_samples = [x.strip() for x in open(input.outliers).readlines()]
        tpm = pandas.read_csv(input.tpm, sep="\t", index_col=0)
        samples = tpm.columns
        tpm_selected = tpm[[sample for sample in samples if sample not in outlier_samples]]
        tpm_selected.index.name = "#"
        sample_data = pandas.read_csv(input.sample_data, sep="\t", index_col=0)
        tpm_selected.columns = [f"ZT{time}" for time in params.timepoints]
        tpm_selected.to_csv(output.tpm, sep="\t")

rule run_ejtk:
    # Sometimes we need to run eJTK first and then run BooteJTK using eJTK's output
    # (For studies with uneven numbers of timepoints
    input:
        "data/{study}/bootejtk/expression.tpm.for_BooteJTK.txt",
    output:
        "data/{study}/ejtk/results.txt",
    resources:
        mem_mb = 6_000,
    shell:
        # Note that the output file name depends upon the number of replicates, which has to vary from one study to the next. Therefore
        # we copy the output to a standardized output name
        "apptainer run --bind {WORKING_DIR} {BOOTEJTK_SIF} /eJTK/eJTK-CalcP.py -f {input} -p /eJTK/ref_files/period24.txt -s /eJTK/ref_files/phases_00-22_by2.txt -a /eJTK/ref_files/asymmetries_02-22_by2.txt -x OUT && cp data/{wildcards.study}/bootejtk/expression.tpm.for_BooteJTK_OUT_jtkout_GammaP.txt {output}"


def bootejtk_rep_count(study):
    timepoints = [x % 24 for x in sample_timepoints(study, drop_outliers=True)]
    counts = list(collections.Counter(timepoints).values())
    if max(counts) != min(counts):
        # BooteJTK gets handled specially if studies have uneven numbers of replicates per timepoints then
        # they also have to simulated as if unique
        return 1
    else:
         # number of replicates per timepoint - same for all timepoints
        return len(timepoints)//len(set(timepoints))

rule run_bootejtk:
    input:
        "data/{study}/expression.tpm.txt",
        "data/{study}/bootejtk/expression.tpm.for_BooteJTK.txt",
        ejtk = lambda wildcards: f"data/{wildcards.study}/ejtk/results.txt" if bootejtk_rep_count(wildcards.study) == 1 else [],
        outliers = lambda wildcards: f"results/{targets[wildcards.study]['tissue']}/outlier_samples.txt",
    output:
        "data/{study}/bootejtk/results.txt",
    resources:
        mem_mb = 6_000,
    params:
        num_reps = lambda wildcards: bootejtk_rep_count(wildcards.study),
        args = lambda wildcards: f'-U -J data/{wildcards.study}/ejtk/results.txt' if bootejtk_rep_count(wildcards.study) == 1 else '',
        out_file_prefix = lambda wildcards: "NoRepSD" if bootejtk_rep_count(wildcards.study) == 1 else 'Vash'
    shell:
        # Note that the output file name depends upon the number of replicates, which has to vary from one study to the next. Therefore
        # we copy the output to a standardized output name
        "apptainer run --bind {WORKING_DIR} {BOOTEJTK_SIF} /BooteJTK/BooteJTK-CalcP.py -f {input[0]} -p /BooteJTK/ref_files/period24.txt -s /BooteJTK/ref_files/phases_00-22_by2.txt -a /BooteJTK/ref_files/asymmetries_02-22_by2.txt -z 25 -r {params.num_reps} -R {params.args} -x OUT && cp data/{wildcards.study}/bootejtk/expression.tpm.for_BooteJTK_{params.out_file_prefix}_OUT_boot25-rep{params.num_reps}_GammaP.txt {output}"

rule plot_qc:
    input:
        salmon_metainfo = lambda wildcards: expand("data/{study}/salmon.meta_info.json", study=studies_by_tissue(wildcards.tissue)),
        expression_tpm = lambda wildcards: expand("data/{study}/label_expression.tpm.txt", study=studies_by_tissue(wildcards.tissue)),
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
    params:
        studies = select_tissue(studies),
    output:
        percent_mapping = "results/{tissue}/qc/percent_mapping.png",
        num_processed = "results/{tissue}/qc/num_processed.png",
        num_mapped = "results/{tissue}/qc/num_mapped.png",
        expression_distributions = "results/{tissue}/qc/expression_distributions.png",
        ENSMUSG00000086503 = "results/{tissue}/qc/Xist_expression.png",
        ENSMUSG00000069036 = "results/{tissue}/qc/Sry_expression.png",
        ENSMUSG00000048376 = "results/{tissue}/qc/F2r_expression.png",
        ENSMUSG00000029368 = "results/{tissue}/qc/Alb_expression.png",
        ENSMUSG00000064337 = "results/{tissue}/qc/mtRnr1_expression.png",
        ENSMUSG00000064339 = "results/{tissue}/qc/mtRnr2_expression.png",
        ENSMUSG00000064336 = "results/{tissue}/qc/mtTf_expression.png",
        ENSMUSG00000054626 = "results/{tissue}/qc/Xlr_expression.png",
        ENSMUSG00000000567 = "results/{tissue}/qc/Sox9_expression.png",
        ENSMUSG00000025332 = "results/{tissue}/qc/Kdm5c_expression.png",
    resources:
        mem_mb = 4000,
    script:
        "scripts/qc_plots.py"

rule plot_genes:
    input:
        expression_tpm = lambda wildcards: expand("data/{study}/label_expression.tpm.txt", study=studies_by_tissue(wildcards.tissue)),
    params:
        studies = select_tissue(studies),
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
        jtk = lambda wildcards: expand("data/{study}/jtk{period}.results.txt",
                                        study=studies_by_tissue(wildcards.tissue),
                                        period=wildcards.period),
        robustness_score = "results/{tissue}/jtk{period}/robustness_score.txt",
    params:
        studies = select_tissue(studies),
    output:
        breakdowns = "results/{tissue}/jtk{period}/breakdowns.png",
        periods = "results/{tissue}/jtk{period}/periods.png",
        amplitudes = "results/{tissue}/jtk{period}/amplitudes.png",
        phases = "results/{tissue}/jtk{period}/phases.png",
        phase_heatmap = "results/{tissue}/jtk{period}/phases.heatmap.loose.png",
        phase_heatmap_robust = "results/{tissue}/jtk{period}/phases.heatmap.robust.png",
        phase_heatmap_svg = "results/{tissue}/jtk{period}/phases.heatmap.loose.svg",
        phase_heatmap_svg_robust = "results/{tissue}/jtk{period}/phases.heatmap.robust.svg",
    script:
        "scripts/plot_jtk.py"

rule plot_overlapped_genes:
    input:
        jtk = lambda wildcards: expand("data/{study}/jtk{period}.results.txt", study=studies_by_tissue(wildcards.tissue), period=wildcards.period),
        jtk_results = "results/{tissue}/jtk{period}.results.txt",
        tpm = lambda wildcards: expand("data/{study}/expression.tpm.txt", study=studies_by_tissue(wildcards.tissue)),
    params:
        studies = select_tissue(studies),
    output:
        num_common_genes = "results/{tissue}/jtk{period}/num_common_genes.txt",
        q_values_overlaps = "results/{tissue}/jtk{period}/q_values_overlaps.txt",
        num_common_genes_heatmap = "results/{tissue}/jtk{period}/num_common_genes_heatmap.png",
        num_qvalue_overlap_heatmap = "results/{tissue}/jtk{period}/num_qvalue_overlap_heatmap.png",
        num_q_p_value_overlap_heatmap = "results/{tissue}/jtk{period}/num_q_pvalue_overlap_heatmap.png",
        common_genes_pvalue = "results/{tissue}/jtk{period}/common_genes_pvalue.txt",
        robustness_score = "results/{tissue}/jtk{period}/robustness_score.txt",
    resources:
        mem_mb = 4000
    script:
        "scripts/plot_overlapped_genes.py"

rule assess_stable_genes:
    input:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
        jtk = "results/{tissue}/jtk24.results.txt",
        spline_fit_summary = "results/{tissue}/spline_fit/summary.full.txt",
    output:
        stable_gene_list = "results/{tissue}/stable_genes/stable_gene_list.txt",
        stable_table = "results/{tissue}/stable_genes/gene_stability.txt",
    script:
        "scripts/assess_stable_genes.py"

rule plot_PCA:
    input:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        jtk = lambda wildcards: expand("data/{study}/jtk.results.txt", study=studies_by_tissue(wildcards.tissue)),
        outlier_samples = "results/{tissue}/outlier_samples.txt",
        robustness = "results/{tissue}/jtk24/robustness_score.txt",
        study_classification = "results/study_classification.json",
        study_table = "results/Liver/study_table.txt"
    params:
        studies = select_tissue(studies),
    resources:
        mem_mb = 6000,
    output:
        dir = directory("results/{tissue}/PCA"),
    script:
        "scripts/plot_PCA.py"

rule plot_robustness:
    input:
        robustness = "results/{tissue}/jtk{period}/robustness_score.txt",
        tpm = lambda wildcards: expand("data/{study}/expression.tpm.txt", study=studies_by_tissue(wildcards.tissue)),
        jtk = lambda wildcards: expand("data/{study}/jtk{period}.results.txt",
                                        study=studies_by_tissue(wildcards.tissue),
                                        period=wildcards.period),
    params:
        studies = select_tissue(studies),
    output:
        expression_level = "results/{tissue}/jtk{period}/robustness/expression_level_robust.png",
        amplitude = "results/{tissue}/jtk{period}/robustness/amplitude_robust.png",
        period = "results/{tissue}/jtk{period}/robustness/period_robust.png",
        phase = "results/{tissue}/jtk{period}/robustness/phase_robust.png",
        histogram = "results/{tissue}/jtk{period}/robustness/histogram.png",
    resources:
        mem_mb = 4000,
    script:
        "scripts/plot_robustness.py"

rule all_samples:
    input:
        tpm = lambda wildcards: expand("data/{study}/expression.tpm.txt", study=studies_by_tissue(wildcards.tissue)),
        num_reads = lambda wildcards: expand("data/{study}/expression.num_reads.txt", study=studies_by_tissue(wildcards.tissue)),
        outliers = lambda wildcards: expand("data/{study}/outlier_samples.txt", study=studies_by_tissue(wildcards.tissue)),
        split_flags = lambda wildcards: expand("data/{study}/samples/split_flag", study=studies_by_tissue(wildcards.tissue)),
    params:
        studies = select_tissue(studies),
        timepoints = lambda wildcards: [sample_timepoints(study) for study in studies_by_tissue(wildcards.tissue)],
    output:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        tpm_parquet = "results/{tissue}/tpm_all_samples.parquet",
        num_reads = "results/{tissue}/num_reads_all_samples.txt",
        num_reads_parquet = "results/{tissue}/num_reads_all_samples.parquet",
        sample_info = "results/{tissue}/all_samples_info.txt",
    resources:
        mem_mb = 6_000,
    run:
        all_tpm = pandas.DataFrame()
        all_num_reads = pandas.DataFrame()
        all_samples = []
        for study, tpmfile in zip(params.studies, input.tpm):
            tpm = pandas.read_csv(tpmfile, sep="\t", index_col="Name")
            time = studies_config.sample_timepoints(study)
            all_samples.append(pandas.DataFrame({
            "study": [study for i in range(len(tpm.columns))],
            "time": time
            }, index=tpm.columns))
            all_tpm = pandas.concat([all_tpm, tpm], axis=1)
        all_tpm.insert(0, 'Symbol', pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol'])
        all_tpm.to_csv(output.tpm, sep ="\t", index="Name")
        all_tpm.to_parquet(output.tpm_parquet)

        all_samples_concat = pandas.concat(all_samples, axis=0)
        all_samples_df = pandas.DataFrame(all_samples_concat)
        all_samples_df.to_csv(output.sample_info, sep="\t")

        for study, numreadsfile in zip(params.studies, input.num_reads):
            num_reads = pandas.read_csv(numreadsfile, sep="\t", index_col="Name")
            all_num_reads = pandas.concat([all_num_reads, num_reads], axis=1)
        all_num_reads.insert(0, 'Symbol', pandas.read_csv("gene_name.txt", sep="\t", index_col="ID")['GeneSymbol'])
        all_num_reads.to_csv(output.num_reads, sep ="\t", index="Name")
        all_num_reads.to_parquet(output.num_reads_parquet)

rule all_jtk: # Gather all JTK results of a tissue together
    input:
        jtk = lambda wildcards: expand("data/{study}/jtk{period}.results.txt",
                                        study=studies_by_tissue(wildcards.tissue),
                                        period=wildcards.period),
    params:
        studies = select_tissue(studies),
    output:
        all_jtk = "results/{tissue}/jtk{period}.results.txt",
    run:
        jtks = []
        for study, jtk_file in zip(params.studies, input.jtk):
            jtk = pandas.read_csv(jtk_file, sep="\t")
            jtk.insert(1, "study", study)
            jtks.append(jtk)
        all_jtk = pandas.concat(jtks, axis=0)
        all_jtk.rename(columns={"CycID": "ID"}, inplace=True)
        all_jtk.to_csv(output.all_jtk, sep="\t", index=False)

checkpoint outlier_detection:
    # NOTE: this wouldn't normally have to be a checkpoint but we want to read its outputs
    # in params and therefore we make it a checkpoint
    input:
        tpm = "data/{study}/expression.tpm.txt",
    output:
        outliers = "data/{study}/outlier_samples.txt",
    script:
        "scripts/outlier_detection.py"

rule gather_outliers:
    input:
        lambda wildcards: expand("data/{study}/outlier_samples.txt", study=studies_by_tissue(wildcards.tissue))
    output:
        "results/{tissue}/outlier_samples.txt"
    shell:
        "cat {input} > {output}"

rule scatter_grid:
    input:
        jtk = lambda wildcards: expand("data/{study}/jtk24.results.txt", study=studies_by_tissue(wildcards.tissue)),
    params:
        studies = select_tissue(studies),
    output:
        amplitude = "results/{tissue}/amplitude_scatter_grid.png",
        phase = "results/{tissue}/phase_scatter_grid.png",
    resources:
        mem_mb = 12000
    script:
        "scripts/plot_scatter_grid.py"

rule plot_consensus_pca:
    input:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
    output:
        outdir = directory("results/{tissue}/consensus_pca/")
    resources:
        mem_mb = 6000
    script:
        "scripts/consensus_pca.py"

rule run_spline_fit:
    input:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
    output:
        summary = "results/{tissue}/spline_fit{permutation}/batches/{batch}.summary.txt",
        curves_fit = "results/{tissue}/spline_fit{permutation}/batches/{batch}.curves_fit.txt",
        curves_pstd = "results/{tissue}/spline_fit{permutation}/batches/{batch}.curves_pstd.txt",
        re = "results/{tissue}/spline_fit{permutation}/batches/{batch}.re.txt",
        re_structure = "results/{tissue}/spline_fit{permutation}/batches/{batch}.re_structure.txt",
        #Note: permutation == '' is  no-permutation
    params:
        num_batches = SPLINE_FIT_N_BATCHES,
    resources:
        mem_mb = 4000
    script:
        "scripts/run_spline_fit.R"

rule gather_spline_fits:
    input:
        summary = expand("results/{{tissue}}/spline_fit{{permutation}}/batches/{batch}.summary.txt",
                            batch=range(SPLINE_FIT_N_BATCHES)),
        curves_fit = expand("results/{{tissue}}/spline_fit{{permutation}}/batches/{batch}.curves_fit.txt",
                            batch=range(SPLINE_FIT_N_BATCHES)),
        curves_pstd = expand("results/{{tissue}}/spline_fit{{permutation}}/batches/{batch}.curves_pstd.txt",
                            batch=range(SPLINE_FIT_N_BATCHES)),
        re = expand("results/{{tissue}}/spline_fit{{permutation}}/batches/{batch}.re.txt",
                            batch=range(SPLINE_FIT_N_BATCHES)),
        re_structure = expand("results/{{tissue}}/spline_fit{{permutation}}/batches/{batch}.re_structure.txt",
                            batch=range(SPLINE_FIT_N_BATCHES)),
    output:
        summary = "results/{tissue}/spline_fit{permutation}/summary.txt",
        curves_fit = "results/{tissue}/spline_fit{permutation}/curves_fit.txt",
        curves_pstd = "results/{tissue}/spline_fit{permutation}/curves_pstd.txt",
        re = "results/{tissue}/spline_fit{permutation}/re.txt",
        re_structure = "results/{tissue}/spline_fit{permutation}/re_structure.txt",
    params:
        num_batches = SPLINE_FIT_N_BATCHES,
        input_dir = "results/{tissue}/spline_fit{permutation}/batches/",
        output_dir = "results/{tissue}/spline_fit{permutation}/",
    resources:
        mem_mb = 6_000
    script:
        "scripts/gather_spline_fits.py"

rule assess_spline_fits:
    input:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
        summary = "results/{tissue}/spline_fit{permutation}/summary.txt",
        curves_fit = "results/{tissue}/spline_fit{permutation}/curves_fit.txt",
        curves_pstd = "results/{tissue}/spline_fit{permutation}/curves_pstd.txt",
        re = "results/{tissue}/spline_fit{permutation}/re.txt",
        re_structure = "results/{tissue}/spline_fit{permutation}/re_structure.txt",
    resources:
        mem_mb = 6_000,
    output:
        out_dir = directory("results/{tissue}/spline_fit{permutation}/stats/"),
        summary = "results/{tissue}/spline_fit{permutation}/summary.full.txt"
    script:
        "scripts/assess_spline_fits.py"

rule plot_spline_fits:
    input:
        tpm = "results/{tissue}/tpm_all_samples.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        summary = "results/{tissue}/spline_fit/summary.full.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
        curves_fit = "results/{tissue}/spline_fit/curves_fit.txt",
        curves_pstd = "results/{tissue}/spline_fit/curves_pstd.txt",
        statsdir = "results/{tissue}/spline_fit/stats/",
    output:
        pca = "results/{tissue}/spline_fit/pca.png",
        tsne = "results/{tissue}/spline_fit/tsne.png",
        phase_distribution = "results/{tissue}/spline_fit/phase_distribution.png",
        phase_correlation = "results/{tissue}/spline_fit/phase_correlation.png",
        gene_plot_dir = directory("results/{tissue}/spline_fit/gene_plots/"),
    resources:
        mem_mb = 4_000,
    script:
        "scripts/plot_spline_fits.py"

rule assess_period_differences:
    input:
        jtk = "results/{tissue}/jtk.results.txt",
        robustness = "results/{tissue}/robustness_score.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
    output:
        "results/{tissue}/assess_jtk/period_statistics.txt",
    script:
        "scripts/assess_period_differences.py"

rule study_table:
    input:
        sample_info = "results/{tissue}/all_samples_info.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
        study_classification = "results/study_classification.json",
    output:
        table = "results/{tissue}/study_table.txt",
    params:
        studies = select_tissue(studies),
    script:
        "scripts/generate_study_table.py"

rule assess_phase_variability:
    input:
        spline_fit_summary = "results/{tissue}/spline_fit/summary.full.txt"
    output:
        phase_std_distribution = "results/{tissue}/spline_fit/phase_variability/phase_std_distribution.png"
    script:
        "scripts/assess_phase_variability.py"

rule run_compare_rhythms:
    input:
        "data/{study1}/expression.num_reads.txt",
        "data/{study2}/expression.num_reads.txt",
        "data/{study1}/sample_data.txt",
        "data/{study1}/expression.tpm.txt",
        "data/{study2}/sample_data.txt",
        "data/{study2}/expression.tpm.txt",
        outliers1 = "data/{study1}/outlier_samples.txt",
        outliers2 = "data/{study2}/outlier_samples.txt",
        split_samples1 = lambda wildcards: checkpoints.split_samples.get(study=wildcards.study1).output[0], # Trigger checkpoint
        split_samples2 = lambda wildcards: checkpoints.split_samples.get(study=wildcards.study2).output[0], # Trigger checkpoint
    output:
        temp("results/{tissue}/compareRhythms/results.{study1}.{study2}.txt")
    params:
        timepoints1 = lambda wildcards: sample_timepoints(wildcards.study1, drop_outliers=True),
        timepoints2 = lambda wildcards: sample_timepoints(wildcards.study2, drop_outliers=True),
    script:
        "scripts/run_compare_rhythms.R"

rule gather_compare_rhythms:
    input:
        compRhythms = lambda wildcards: [f"results/{wildcards.tissue}/compareRhythms/results.{study1}.{study2}.txt"
                                            for study1 in studies_by_tissue(wildcards.tissue)
                                            for study2 in studies_by_tissue(wildcards.tissue)
                                            if study1 > study2]
    output:
        "results/{tissue}/compareRhythms/summary.txt"
    run:
        import pandas
        summary = []
        for study1 in studies_by_tissue(wildcards.tissue):
            for study2 in studies_by_tissue(wildcards.tissue):
                if study1 <= study2:
                    continue
                res = pandas.read_csv(f"results/{wildcards.tissue}/compareRhythms/results.{study1}.{study2}.txt", sep="\t")
                table = pandas.DataFrame(dict(
                    counts = res.category.value_counts(),
                    study1 = study1,
                    study2 = study2,
                ))
                table.index.name = "category"
                summary.append(table)
        pandas.concat(summary).to_csv(output[0], sep="\t")

rule plot_compare_rhythms:
    input:
        summary = "results/{tissue}/compareRhythms/summary.txt",
        all_jtk = "results/{tissue}/jtk24.results.txt",
    output:
        outdir = directory("results/{tissue}/compareRhythms/plots/")
    script:
        "scripts/plot_compare_rhythms.py"

rule prepare_supplemental:
    input:
        jtk = "results/{tissue}/jtk24.results.txt",
        robustness_score = "results/{tissue}/jtk24/robustness_score.txt",
        sample_info = "results/{tissue}/all_samples_info.txt",
        study_table = "results/{tissue}/study_table.txt",
        tpm = "results/{tissue}/tpm_all_samples.txt",
        num_reads = "results/{tissue}/num_reads_all_samples.txt",
        outliers = "results/{tissue}/outlier_samples.txt",
        sim_summary = "results/{tissue}/spline_fit/summary.full.txt",
        sim_curves_fit = "results/{tissue}/spline_fit/curves_fit.txt",
        sim_curves_pstd = "results/{tissue}/spline_fit/curves_pstd.txt",
        sim_re = "results/{tissue}/spline_fit/re.txt",
        sim_re_structure = "results/{tissue}/spline_fit/re_structure.txt",
    output:
        outdir = directory("results/{tissue}/supplemental/")
    script:
        "scripts/prepare_supplemental.py"
