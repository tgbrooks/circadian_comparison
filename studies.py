import re
import pandas
def extract_ctzt(times):
    return [int(re.search("[ZC]T(\d+)", time).groups()[0]) for time in times]
def manella21_time(times):
    return [int(re.search("(\d+)[AB]", time).groups()[0]) for time in times]
def only_number(times):
    return [int(re.search("(\d+)", time).groups()[0]) for time in times]

def sample_timepoints(study):
    sample_data = pandas.read_csv(f"data/{study}/sample_data.txt", sep="\t", index_col="geo_accession")
    expression_table = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0, nrows=5)
    times = targets[study]["time"](sample_data, expression_table)
    return times

targets = {
    "Schwartz21": {
        "GSE": "GSE165198",
        "sample_selector": lambda x: True,
    },

    "Yang16A_M": {
        "GSE": "GSE70497",
        "sample_selector": lambda x: x['genotype/variation'] == "WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['sample collection time'])),
        "tissue": "Liver",
    },

    "Yang16B_M": {
        "GSE": "GSE70499",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.Sex == "male",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['sample collection time'])),
        "tissue": "Liver",
    },

    "Yang16B_F": {
        "GSE": "GSE70499",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.Sex == "female",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['sample collection time'])),
        "tissue": "Liver",
    },

    "Lahens15": {
        "GSE": "GSE40190",
        "sample_selector": lambda x: True,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].time)),
        "tissue": "Liver",
    },

    "Weger18": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
    },

    "Weger18_Liver_M": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "male" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
    },

    "Weger18_Liver_F": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "Female" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
     },

    "Zhang14_RNAseq_Liver_M": {
        "GSE": "GSE54651",
        "sample_selector": lambda x: x.tissue == "liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
    },

    "Pan19": {
        "GSE": "GSE130890",
        "sample_selector": lambda x: x["genotype/variation"] =="WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
    },

    "Morton20_Liver": {
        "GSE": "GSE151565",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.tissue == "Liver",
        "time": lambda sample_data, expression_table: only_number(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
    },

    "Atger15_AdLib": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Ad Libitum",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
    },

    "Atger15_NightFeed": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Night restricted",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
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
        "tissue": "Liver",
    },

    "Guan20_Liver": {
        "GSE": "GSE143524",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.tissue == "Liver" and x.description == "Ad_lib_feeding",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
    },

    "Koronowski19_F": {
        "GSE": "GSE117134",
        "sample_selector": lambda x: x['genotype/variation'] == "WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time of harvest'])),
        "tissue": "Liver",
    },

    "Meng20": {
        "GSE": "GSE150888",
        "sample_selector": lambda x: x.genotype == "Xbp1 flx/flx (WT, wild-type)",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
    },

    "Xin21_Liver_NightFeed": {
        "GSE": "GSE150380",
        "sample_selector": lambda x: x['dietary regimen'] == "NRF",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time in hours'])),
        "tissue": "Liver",
    },

    "Yang20": {
        "GSE": "GSE115264",
        "sample_selector": lambda x: x.genotype == "WT",
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['sample collection time']),
        "tissue": "Liver",
    },

    "Kinouchi18_Liver": {
        "GSE": "GSE107787",
        "sample_selector": lambda x: x.source_name_ch1 == "Liver" and "ad libitum" in x['timepoint/condition'],
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
    },

    "Weger20_Bmal1WT": {
        "GSE": "GSE135898",
        "sample_selector": lambda x: x.genotype == "Bmal1  WT" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
    },

    "Weger20_HlfDbpTefWT": {
        "GSE": "GSE135875",
        "sample_selector": lambda x: x.genotype == "Hlf+/+/Dbp+/+/Tef+/+" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
    },

    "Weger20_CryWT": {
        "GSE": "GSE135898",
        "sample_selector": lambda x: x.genotype == "Cry1/2_WT" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
    },

    "Mermet18": {
        "GSE": "GSE101423",
        #"GSE_postfix": "-GPL17021",
        "GSE_postfix": "-GPL19057",
        "sample_selector":  lambda x: x.genotype == "Cry1IntronWT" and x.tissue == "Liver" and x.source_name_ch1 == "PolyA-selected RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },
}

# List of studies to perform
studies = ["Lahens15", "Weger18_Liver_M", "Weger18_Liver_F", "Zhang14_RNAseq_Liver_M", "Pan19", "Morton20_Liver", "Atger15_AdLib", "Atger15_NightFeed", "Manella21_Liver", "Weger18", "Yang16A_M", "Yang16B_M", "Yang16B_F", "Guan20_Liver", "Koronowski19_F", "Xin21_Liver_NightFeed", "Yang20", "Kinouchi18_Liver", "Weger20_Bmal1WT", "Weger20_HlfDbpTefWT", "Weger20_CryWT", "Mermet18"]

# List of all tissues being run
tissues = list(set(targets[study]['tissue'] for study in studies))
