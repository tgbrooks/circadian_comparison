import re
import pandas
def extract_ctzt(times):
    return [int(re.search("[ZC]T(\d+)", time).groups()[0]) for time in times]
def manella21_time(times):
    return [int(re.search("(\d+)[AB]", time).groups()[0]) for time in times]
def only_number(times):
    return [int(re.search("(\d+)", time).groups()[0]) for time in times]
def hirako18_time(times):
    # Times are in clock time with ZT0 = 7:00am, so 0 Hour (midnight) is at ZT17
    return [int(re.search("(\d+)", time).groups()[0])+17 for time in times]

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
        "highlight": True,
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

    "Weger18_Liver_M": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "male" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "highlight": True,
    },

    "Weger18_Liver_F": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "Female" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
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
        "highlight": True,
    },

    "Atger15_AdLib": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Ad Libitum",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "highlight": True,
    },

    "Atger15_NightFeed": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Night restricted",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq":  "RiboZero",
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
        "seq": "3prime",
        "highlight": True,
    },

    "Guan20_Liver": {
        "GSE": "GSE143524",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.tissue == "Liver" and x.description == "Ad_lib_feeding",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "highlight": True,
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
        "seq":  "RiboZero",
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
        "highlight": True,
    },

    "Weger20_HlfDbpTefWT": {
        "GSE": "GSE135875",
        "sample_selector": lambda x: x.genotype == "Hlf+/+/Dbp+/+/Tef+/+" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Weger20_CryWT": {
        "GSE": "GSE135898",
        "sample_selector": lambda x: x.genotype == "Cry1/2_WT" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Mermet18": {
        "GSE": "GSE101423",
        #"GSE_postfix": "-GPL17021",
        "GSE_postfix": "-GPL19057",
        "sample_selector":  lambda x: x.genotype == "Cry1IntronWT" and x.tissue == "Liver" and x.source_name_ch1 == "PolyA-selected RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },
#start
    "Benegiamo18": {
        "GSE": "GSE98042",
        "sample_selector": lambda x: x.genotype == "wild type" and x.tissue == "total liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
    },

    "Cajan16": {
        "GSE": "GSE61775",
        "GSE_postfix": "-GPL17021",
        "sample_selector": lambda x: "Liver" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },

    "Chaix19_AdLib_HFD": {
        "GSE": "GSE102072",
        "GSE_postfix": "-GPL19057", # note: only ad lib feeding WTs are in this series matrix
        "sample_selector": lambda x: x.genotype == "WT" and x['feeding group'] == "FA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
    },

    "Chaix19_NightFeed_HFD": {
        "GSE": "GSE102072",
        "GSE_postfix": "-GPL17021", # note: only TRF WTs are in this series matrix
        "sample_selector": lambda x: x.genotype == "WT" and x['feeding group'] == "FT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
    },

    "Chen19": {
        "GSE": "GSE133342",
        "sample_selector": lambda x: "CT" in x.title,
        "time": lambda sample_data, expression_table: only_number(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },

    "Du14": {
        "GSE": "GSE57313",
        "sample_selector": lambda x: x.genotype == "control littermate",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
    },

    "Fader19": {
        "GSE": "GSE119780",
        "sample_selector": lambda x: x.treatment == "Sesame Oil Control",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Gaucher19_Chronic_Cntrl": {
        "GSE": "GSE132103",
        "sample_selector": lambda x: "CHRONIC_CTRL" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },

    "Gaucher19_Acute_Cntrl": {
        "GSE": "GSE132103",
        "sample_selector": lambda x: "ACUTE_CTRL" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        # NOTE: this study has irregular timepoint CT1 resulting in run JTK error
    },

    "Greenwell19_AdLib": {
        "GSE": "GSE118967",
        "sample_selector": lambda x: "LB" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "3prime",
        "highlight": True,
    },

    "Greenwell19_NightFeed": {
        "GSE": "GSE118967",
        "sample_selector": lambda x: "NR" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "3prime",
    },

    "Hidalgo19": {
        "GSE": "GSE125867",
        "sample_selector": lambda x: x.genotype == "WT",
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['zeitgeber time']),
        "tissue": "Liver",
    },

    "Hirako18": {
        "GSE": "GSE109908",
        "sample_selector": lambda x: x.agent == "Control",
        "time": lambda sample_data, expression_table: hirako18_time(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
    },

    "Janich15": {
        "GSE": "GSE67305",
        "sample_selector": lambda x: x['insert type'] == "total RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "tissue": "Liver",
    },

    "Levine20": {
        "GSE": "GSE133989",
        "sample_selector": lambda x: x.genotype == "WT" and x.treatment == "H2O" and x.library_source == "transcriptomic",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['collection time point'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Li19_Young": {
        "GSE": "GSE113745",
        "sample_selector": lambda x: "Young" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
    },

    "Li19_Old": {
        "GSE": "GSE113745",
        "sample_selector": lambda x: "Old" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
    },

    "Menet12": {
        "GSE": "GSE36871",
        "sample_selector": lambda x: x.tissue == "Liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Quagliarini19_NormalDiet": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "control REP" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Quagliarini19_HFD": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "HFD REP" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },

    "Quagliarini19_NormalDiet_WTvsKO": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "LFD WT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "highlight": True,
    },

    "Quagliarini19_HFD_WTvsKO": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "HFD WT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
    },

    "Stubblefield18": {
        "GSE": "GSE105413",
        "sample_selector": lambda x: x.diet == "Regular Chow" and "ZT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "highlight": True,
    },

    "Wu19": {
        "GSE": "GSE138019",
        "sample_selector": lambda x: x.genotype == "WT" and "ZT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
    },

    "Yeung17": {
        "GSE": "GSE100457",
        "GSE_postfix": "-GPL17021",
        "sample_selector": lambda x: x.genotype == "WT" and "RNASeq" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Kidney",
    },

    "Zhang14_RNAseq_Kidney_M": {
        "GSE": "GSE54651",
        "sample_selector": lambda x: x.tissue == "kidney",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Kidney",
    },

    "Castelo-Szekely17": {
        "GSE": "GSE81283",
        "sample_selector": lambda x: "total RNA" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Kidney",
    },
}

# List of studies to perform
studies = [
    # Highlighted studies -- all consistent design
    "Weger18_Liver_M", "Morton20_Liver", "Atger15_AdLib", "Manella21_Liver", "Yang16B_M", "Weger20_Bmal1WT", "Weger20_HlfDbpTefWT", "Weger20_CryWT", "Fader19", "Greenwell19_AdLib", "Levine20", "Quagliarini19_NormalDiet", "Quagliarini19_NormalDiet_WTvsKO", "Stubblefield18", "Menet12", "Guan20_Liver",
    # Possible highlights
    "Hidalgo19", "Li19_Young", "Li19_Old", "Du14", "Pan19",
    # Remaining
    "Lahens15", "Weger18_Liver_F", "Zhang14_RNAseq_Liver_M", "Atger15_NightFeed", "Yang16A_M", "Yang16B_F", "Koronowski19_F", "Xin21_Liver_NightFeed", "Yang20", "Kinouchi18_Liver", "Mermet18", "Benegiamo18", "Cajan16", "Chaix19_AdLib_HFD", "Chaix19_NightFeed_HFD", "Chen19", "Gaucher19_Chronic_Cntrl", "Greenwell19_NightFeed", "Hirako18", "Quagliarini19_HFD", "Quagliarini19_HFD_WTvsKO", "Wu19", "Yeung17", "Zhang14_RNAseq_Kidney_M", "Castelo-Szekely17", ]
#studies = ["Lahens15", "Weger18_Liver_M", "Weger18_Liver_F", "Zhang14_RNAseq_Liver_M", "Pan19", "Morton20_Liver", "Atger15_AdLib", "Atger15_NightFeed", "Manella21_Liver", "Yang16A_M", "Yang16B_M", "Yang16B_F", "Guan20_Liver", "Koronowski19_F", "Xin21_Liver_NightFeed", "Yang20", "Kinouchi18_Liver", "Weger20_Bmal1WT", "Weger20_HlfDbpTefWT", "Weger20_CryWT", "Mermet18", "Benegiamo18", "Cajan16", "Chaix19_AdLib_HFD", "Chaix19_NightFeed_HFD", "Chen19", "Du14", "Fader19", "Gaucher19_Chronic_Cntrl", "Greenwell19_AdLib", "Greenwell19_NightFeed", "Hidalgo19", "Hirako18", "Levine20", "Li19_Young", "Li19_Old", "Menet12", "Quagliarini19_NormalDiet", "Quagliarini19_HFD", "Quagliarini19_NormalDiet_WTvsKO", "Quagliarini19_HFD_WTvsKO", "Stubblefield18", "Wu19", "Yeung17", "Zhang14_RNAseq_Kidney_M", "Castelo-Szekely17", ]

# List of all tissues being run
tissues = list(set(targets[study]['tissue'] for study in studies))

def studies_by_tissue(tissue):
    return [study for study in studies if targets[study]['tissue'] == tissue]

def select_tissue(studies):
    def f(wildcards, output):
        return [study for study in studies if targets[study]['tissue'] == wildcards.tissue]
    return f
