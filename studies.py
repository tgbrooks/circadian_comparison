import re
import pandas
def extract_ctzt(times):
    return [int(re.search("[ZC]T(\d+)", time).groups()[0]) for time in times]
def manella21_time(times):
    return [int(re.search("(\d+)[AB]", time).groups()[0]) for time in times]
def morton20_time(times):
    return [int(re.search("(\d+)", time).groups()[0]) for time in times]

def sample_timepoints(study):
    sample_data = pandas.read_csv(f"data/{study}/sample_data.txt", sep="\t", index_col="geo_accession")
    expression_table = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0)
    times = targets[study]["time"](sample_data, expression_table)
    return times

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
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].time)),
    },

    "Weger18": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
    },

    "Weger18_Liver_M": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "male" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
    },

    "Weger18_Liver_F": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "Female" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
     },

    "Zhang14_RNAseq_Liver_M": {
        "GSE": "GSE54651",
        "sample_selector": lambda x: x.tissue == "liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
    },

    "Pan19": {
        "GSE": "GSE130890",
        "sample_selector": lambda x: x["genotype/variation"] =="WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
    },

    "Morton20_Liver": {
        "GSE": "GSE151565",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.tissue == "Liver",
        "time": lambda sample_data, expression_table: morton20_time(list(sample_data.loc[expression_table.columns]['time point'])),
    },

    "Atger15_AdLib": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Ad Libitum",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
    },

    "Atger15_NightFeed": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Night restricted",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
    },

    "Koike12_RNAseq": {
        "GSE": "GSE39978",
        "sample_selector": lambda x: x["genotype/variation"] == "wild-type",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].time)),
        # NOTE: this study has SOLiD data, not compatible with Illumina data. Would need another way to align
    },

    "Manella21_Liver": {
        "GSE": "GSE159135",
        "sample_selector": lambda x: x.tissue == "Liver" and x.genotype == "Alb-Cre" and x["feeding regimen"] == "AL",
        "time": lambda sample_data, expression_table: manella21_time(list(sample_data.loc[expression_table.columns].title)),
    },

}

# List of studies to perform
studies = ["Lahens15", "Weger18_Liver_M", "Weger18_Liver_F", "Zhang14_RNAseq_Liver_M", "Pan19", "Morton20_Liver", "Atger15_AdLib", "Atger15_NightFeed", "Manella21_Liver"]
