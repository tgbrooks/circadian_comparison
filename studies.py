import re
import pandas
def extract_ctzt(times):
    return [int(re.search("[ZC]T(\d+)", time).groups()[0]) for time in times]
def extract_ctzt_lower(times):
    return [int(re.search("[zc]t(\d+)", time).groups()[0]) for time in times]
def extract_gap_ctzt(times):
    return [int(re.search("[ZC]T (\d+)", time).groups()[0]) for time in times]
def extract_TP(times):
    return [int(re.search("TP(\d+)", time).groups()[0]) for time in times]
def manella21_time(times):
    return [int(re.search("(\d+)[AB]", time).groups()[0]) for time in times]
def only_number(times):
    return [int(re.search("(\d+)", time).groups()[0]) for time in times]
def hirako18_time(times):
    # Times are in clock time with ZT0 = 7:00am, so 0 Hour (midnight) is at ZT17
    return [int(re.search("(\d+)", time).groups()[0])+17 for time in times]

def sample_timepoints(study, drop_outliers=False):
    sample_data = pandas.read_csv(f"data/{study}/sample_data.txt", sep="\t", index_col="geo_accession")
    expression_table = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0, nrows=5)
    times = targets[study]["time"](sample_data, expression_table)
    if drop_outliers:
        outlier_samples = [x.strip() for x in open(f"results/{targets[study]['tissue']}/outlier_samples.txt").readlines()]
        return list(time for time, sample in zip(times, expression_table.columns) if sample not in outlier_samples)
    return times

def study_samples(study):
    expression_table = pandas.read_csv(f"data/{study}/expression.tpm.txt", sep="\t", index_col=0, nrows=5)
    return expression_table.columns

def study_samples_is_outliers(stuy, tissue):
    samples = study_samples(study)
    outliers = [x.strip() for x in open(f"results/{tissue}/outlier_samples.txt", "r").read().strip().splitlines()]
    return [(s in outliers) for s in samples]

targets = {
    #"Schwartz21": {
    #    "GSE": "GSE165198",
    #    "sample_selector": lambda x: True,
    #},

    #### LIVER ####
    "Yang16A_M": {
        "GSE": "GSE70497",
        "sample_selector": lambda x: x['genotype/variation'] == "WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['sample collection time'])),
        "tissue": "Liver",
        "short_name": "Yang16A",
        "seq": "PolyA",
        "sex": "M",
        "age_low": 16,
        "age_high": 24,
        "light": "DD",
        "D_time": 36,
        "note": "Bmal1 f/f; tamoxifen treated",
        "PMID": ['26843191'],
    },

    "Yang16B_M": {
        "GSE": "GSE70499",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.Sex == "male",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['sample collection time'])),
        "tissue": "Liver",
        "highlight": True,
        "seq": "PolyA",
        "short_name": "Yang16B",
        "sex": "M",
        "age_low": 6,
        "age_high": 14,
        "light": "LD",
        "PMID": ['26843191'],
    },

    "Yang16B_F": {
        "GSE": "GSE70499",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.Sex == "female",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['sample collection time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Yang16C",
        "sex": "F",
        "age_low": 6,
        "age_high": 14,
        "light": "LD",
        "PMID": ['26843191'],
    },

    "Lahens15": {
        "GSE": "GSE40190",
        "sample_selector": lambda x: True,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].time)),
        "tissue": "Liver",
        "seq": "PolyA",
        "sex": "M",
        "light": "DD",
        "D_time": 36,
        "age_low": 6,
        "age_high": 6,
        "PMID": [],
    },

    "Weger18_Liver_M": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "male" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "highlight": True,
        "short_name": "Weger19A",
        "sex": "M",
        "light": "LD",
        "age_low": 15,
        "age_high": 16,
    },

    "Weger18_Liver_F": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "Female" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "short_name": "Weger19B",
        "sex": "F",
        "light": "LD",
        "age_low": 15,
        "age_high": 16,
     },

    "Zhang14_RNAseq_Liver_M": {
        "GSE": "GSE54651",
        "sample_selector": lambda x: x.tissue == "liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Zhang14",
        "sex": "M",
        "light": "DD",
        "D_time": 30,
        "age_low": 6,
        "age_high": 7,
    },

    "Pan19": {
        "GSE": "GSE130890",
        "sample_selector": lambda x: x["genotype/variation"] =="WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "seq": "PolyA",
        "tissue": "Liver",
        "sex": "M",
        "short_name": "Pan20",
        "light": "DD",
        "D_time": 24,
        "age_low": 8,
        "age_high": 12,
        "note": "XBP1 flox",
    },

    "Morton20_Liver": {
        "GSE": "GSE151565",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.tissue == "Liver",
        "time": lambda sample_data, expression_table: only_number(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Morton20",
        "sex": "M",
        "light": "LD",
        "age_low": 26,
        "age_high": 26,
    },

    "Atger15_AdLib": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Ad Libitum",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "highlight": True,
        "short_name": "Sinturel17A",
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 14,
    },

    "Atger15_NightFeed": {
        "GSE": "GSE73552",
        "sample_selector": lambda x: x.genotype =="Wild type" and x.feeding == "Night restricted",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "short_name": "Sinturel17B",
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 14,
        "note": "night-restricted feeding"
    },

    #"Koike12_RNAseq": {
    #    "GSE": "GSE39978",
    #    "sample_selector": lambda x: x["genotype/variation"] == "wild-type",
    #    "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].time)),
    #    # NOTE: this study has SOLiD data, not compatible with Illumina data. Would need another way to align
    #    "short_name": "Koike12",
    #    "seq": "SOLiD",
    #},

    "Manella21_Liver": {
        "GSE": "GSE159135",
        "sample_selector": lambda x: x.tissue == "Liver" and x.genotype == "Alb-Cre" and x["feeding regimen"] == "AL",
        "time": lambda sample_data, expression_table: manella21_time(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq": "3prime",
        "short_name": "Manella21",
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 16,
        "note": "Alb-Cre+",
        #"highlight": True, # Genotype is not strictly WT
    },

    "Guan20_Liver": {
        "GSE": "GSE143524",
        "sample_selector": lambda x: x['genotype/variation'] == "WT" and x.tissue == "Liver" and x.description == "Ad_lib_feeding",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Guan20",
        "sex": "M",
        "light": "LD",
        "age_low": 8,
        "age_high": 12,
        "note": "Rev-erb alpha fl/fl Rev-erb beta fl/fl",
    },

    "Koronowski19_F": {
        "GSE": "GSE117134",
        "sample_selector": lambda x: x['genotype/variation'] == "WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time of harvest'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "short_name": "Mortimer21",
        "sex": "F",
        "light": "LD",
        "age_low": 8,
        "age_high": 12,
        "note": "Alfp-Cre -/tg",
    },

    "Meng20": {
        "GSE": "GSE150888",
        "sample_selector": lambda x: x.genotype == "Xbp1 flx/flx (WT, wild-type)",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "sex": "M",
        "light": "DD",
        "age_low": 12,
        "age_high": 16,
    },

    "Xin21_Liver_NightFeed": {
        "GSE": "GSE150380",
        "sample_selector": lambda x: x['dietary regimen'] == "NRF",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time in hours'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "short_name": "Xin21",
        "sex": "F",
        "light": "LD",
        "age_low": 9,
        "age_high": 9,
        "note": "night-restricted feeding"
    },

    "Yang20": {
        "GSE": "GSE115264",
        "sample_selector": lambda x: x.genotype == "WT",
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['sample collection time']),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Brooks22",
        "sex": "M",
        "light": "DD",
        "age_low": 16,
        "age_high": 24,
        "D_time": 36,
        "note": "Bmal1 fx/fx",
    },

    "Kinouchi18_Liver": {
        "GSE": "GSE107787",
        "sample_selector": lambda x: x.source_name_ch1 == "Liver" and "ad libitum" in x['timepoint/condition'],
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Kinouchi18",
        "sex": "M",
        "light": "LD",
        "age_low": 8,
        "age_high": 8,
    },

    "Weger20_Bmal1WT": {
        "GSE": "GSE135898",
        "sample_selector": lambda x: x.genotype == "Bmal1  WT" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Weger21A",
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 14,
    },

    "Weger20_HlfDbpTefWT": {
        "GSE": "GSE135875",
        "sample_selector": lambda x: x.genotype == "Hlf+/+/Dbp+/+/Tef+/+" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Weger21B",
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 14,
    },

    "Weger20_CryWT": {
        "GSE": "GSE135898",
        "sample_selector": lambda x: x.genotype == "Cry1/2_WT" and x['feeding regimen'] == "AL",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Weger21C",
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 14,
    },

    "Mermet18": {
        "GSE": "GSE101423",
        #"GSE_postfix": "-GPL17021",
        "GSE_postfix": "-GPL19057",
        "sample_selector":  lambda x: x.genotype == "Cry1IntronWT" and x.tissue == "Liver" and x.source_name_ch1 == "PolyA-selected RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "sex": "M",
        "light": "LD",
        "age_low": 8,
        "age_high": 12,
    },

    "Benegiamo18": {
        "GSE": "GSE98042",
        "sample_selector": lambda x: x.genotype == "wild type" and x.tissue == "total liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "sex": "M",
        "light": "LD",
        "age_low": 14,
        "age_high": 14,
    },

    "Cajan16": {
        "GSE": "GSE61775",
        "GSE_postfix": "-GPL17021",
        "sample_selector": lambda x: "Liver" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "sex": "M",
        "light": "LD",
        "age_low": 10,
        "age_high": 12,
    },

    "Chaix19_AdLib_HFD": {
        "GSE": "GSE102072",
        "GSE_postfix": "-GPL19057", # note: only ad lib feeding WTs are in this series matrix
        "sample_selector": lambda x: x.genotype == "WT" and x['feeding group'] == "FA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Chaix19A",
        "sex": "M",
        "light": "LD",
        "age_low": 24,
        "age_high": 24,
        "note": "high fat diet",
    },

    "Chaix19_NightFeed_HFD": {
        "GSE": "GSE102072",
        "GSE_postfix": "-GPL17021", # note: only TRF WTs are in this series matrix
        "sample_selector": lambda x: x.genotype == "WT" and x['feeding group'] == "FT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Chaix19B",
        "sex": "M",
        "light": "LD",
        "age_low": 24,
        "age_high": 24,
        "note": "night-restricted feeding; high fat diet"
    },

    "Chen19": {
        "GSE": "GSE133342",
        "sample_selector": lambda x: "CT" in x.title,
        "time": lambda sample_data, expression_table: only_number(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "short_name": "Li20",
        "sex": "M",
        "light": "DD",
        "D_time": 1008,
        "age_low": 6,
        "age_high": 6,
        "note": "six weeks of DD",
    },

    "Du14": {
        "GSE": "GSE57313",
        "sample_selector": lambda x: x.genotype == "control littermate",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 24,
        "note": "Alb Cr-ERT2; tamoxifen treated",
    },

    "Fader19": {
        "GSE": "GSE119780",
        "sample_selector": lambda x: x.treatment == "Sesame Oil Control",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 10,
        "note": "sesame oil gavage",
    },

    "Gaucher19_Chronic_Cntrl": {
        "GSE": "GSE132103",
        "sample_selector": lambda x: "CHRONIC_CTRL" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Gaucher19",
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 24,

    },

    "Gaucher19_Acute_Cntrl": {
        "GSE": "GSE132103",
        "sample_selector": lambda x: "ACUTE_CTRL" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Gaucher19B",
        # NOTE: this study has irregular timepoint CT1 resulting in run JTK error
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 24,
    },

    "Greenwell19_AdLib": {
        "GSE": "GSE118967",
        "sample_selector": lambda x: "LB" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "3prime",
        "highlight": True,
        "short_name": "Greenwell19A",
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 13,
    },

    "Greenwell19_NightFeed": {
        "GSE": "GSE118967",
        "sample_selector": lambda x: "NR" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "3prime",
        "short_name": "Greenwell19B",
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 13,
        "note": "night-restricted feeding"
    },

    "Hidalgo19": {
        "GSE": "GSE125867",
        "sample_selector": lambda x: x.genotype == "WT",
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['zeitgeber time']),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Rubio-Ponce21",
        "highlight": True,
        "sex": "M",
        "light": "LD",
        "age_low": 8,
        "age_high": 9,
    },

    "Hirako18": {
        "GSE": "GSE109908",
        "sample_selector": lambda x: x.agent == "Control",
        "time": lambda sample_data, expression_table: hirako18_time(list(sample_data.loc[expression_table.columns]['time point'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "sex": "F",
        "light": "LD",
        "age_low": 10,
        "age_high": 10,
    },

    "Janich15": {
        "GSE": "GSE67305",
        "sample_selector": lambda x: x['insert type'] == "total RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "tissue": "Liver",
        "seq": "unknown",
        "trim_adaptor": "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --match-read-wildcards -m 6",
        "sex": "M",
        "light": "LD",
        "age_low": 11,
        "age_high": 12,
    },

    "Levine20": {
        "GSE": "GSE133989",
        "sample_selector": lambda x: x.genotype == "WT" and x.treatment == "H2O" and x.library_source == "transcriptomic",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['collection time point'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "sex": "M",
        "light": "LD",
        "age_low": 32,
        "age_high": 32,
    },

    "Li19_Young": {
        "GSE": "GSE113745",
        "sample_selector": lambda x: "Young" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Li19A",
        "sex": "M",
        "light": "LD",
        "age_low": 16,
        "age_high": 16,
    },

    "Li19_Old": {
        "GSE": "GSE113745",
        "sample_selector": lambda x: "Old" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Li19B",
        "sex": "M",
        "light": "LD",
        "age_low": 76,
        "age_high": 76,
    },

    "Menet12": {
        "GSE": "GSE36871",
        "sample_selector": lambda x: x.tissue == "Liver",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Trott18",
        "highlight": True,
        "sex": "M",
        "light": "LD",
        "age_low": 12,
        "age_high": 24,
    },

    "Quagliarini19_NormalDiet": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "control REP" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Quagliarini19A",
        "sex": "M",
        "light": "LD",
        "age_low": 17,
        "age_high": 18
    },

    "Quagliarini19_HFD": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "HFD REP" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Quagliarini19B",
        "sex": "M",
        "light": "LD",
        "age_low": 17,
        "age_high": 18,
        "note": "high fat diet",
    },

    "Quagliarini19_NormalDiet_WTvsKO": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "LFD WT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "highlight": True,
        "short_name": "Quagliarini19C",
        "sex": "M",
        "light": "LD",
        "age_low": 17,
        "age_high": 18,
    },

    "Quagliarini19_HFD_WTvsKO": {
        "GSE": "GSE108688",
        "sample_selector": lambda x: "HFD WT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Liver",
        "seq": "PolyA",
        "short_name": "Quagliarini19D",
        "sex": "M",
        "light": "LD",
        "age_low": 17,
        "age_high": 18,
        "note": "high fat diet",
    },

    "Stubblefield18": {
        "GSE": "GSE105413",
        "sample_selector": lambda x: x.diet == "Regular Chow" and "ZT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "highlight": True,
        "sex": "M",
        "light": "LD",
        "age_low": 9,
        "age_high": 12,
    },

    "Wu19": {
        "GSE": "GSE138019",
        "sample_selector": lambda x: x.genotype == "WT" and "ZT" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber'])),
        "tissue": "Liver",
        "seq":  "RiboZero",
        "sex": "M",
        "light": "LD",
        "age_low": 16,
        "age_high": 24,
    },

    "Aviram21": {
        "GSE": "GSE171975",
        "tissue": "Liver",
        "sample_selector": lambda x: x.genotype == "WT",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "seq": "3prime",
        "sex": "M",
        "light": "DD",
        "age_low": 12,
        "age_high": 12,
    },

    "Frazier22A": {
        "GSE": "GSE184303",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.genotype == 'Bmal1 fl/fl') and (x['microbial status'] == 'SPF'),
        "time": lambda sample_data, expression_table: extract_gap_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "seq": "PolyA",
        "sex": "M",
        "light": "LD",
        "age_low": 13,
        "age_high": 17,
        "note": "Bmal1 fl/fl; SPF", # Specific Pathogen Free
    },

    "Frazier22B": {
        "GSE": "GSE184303",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.genotype == 'Bmal1 fl/fl') and (x['microbial status'] == 'GF'),
        "time": lambda sample_data, expression_table: extract_gap_ctzt(list(sample_data.loc[expression_table.columns]['timepoint'])),
        "seq": "PolyA",
        "sex": "M",
        "light": "LD",
        "age_low": 13,
        "age_high": 17,
        "note": "Bmal1 fl/fl; GF", # Germ free
    },

    "Wu23": {
        "GSE": "GSE195456",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.genotype == 'wild type') and (x.time not in ['ZT(Hightest)', 'ZT(Lowest)']),
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "seq": "PolyA",
        "sex": "M",
        "light": "LD",
        "age_low": 8,
        "age_high": 52,
    },

    "Acosta-Rodriguez22A":  {
        "GSE": "GSE190939",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.age == '6 months-old') and (x.feeding == 'Ad libitum'),
        "time": lambda sample_data, expression_table: manella21_time(list(sample_data.loc[expression_table.columns]['title'])),
        "short_name": "Acosta-Rodríguez22A",
        "seq": "PolyA",
        "sex": "M",
        "ligh": "DD",
        "age_low": 26,
        "age_high": 26,
    },

    "Acosta-Rodriguez22B": {
        "GSE": "GSE190939",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.age == '19 months-old') and (x.feeding == 'Ad libitum'),
        "time": lambda sample_data, expression_table: manella21_time(list(sample_data.loc[expression_table.columns]['title'])),
        "short_name": "Acosta-Rodríguez22B",
        "seq": "PolyA",
        "sex": "M",
        "light": "DD",
        "age_low": 82,
        "age_high": 82,
    },

    "Katsioudi22": {
        "GSE": "GSE208768",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.genotype == 'WT'),
        "time": lambda sample_data, expression_table: only_number(list(sample_data.loc[expression_table.columns]['zeitgeber time'])),
        "seq": "RiboZero",
        "sex": "M",
        "light": "LD",
        "age_low": 52,
        "age_high": 104,
    },

    "Astafev23": {
        "GSE": "GSE216416",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.treatment == 'Ad-libitum diet'),
        "time": lambda sample_data, expression_table: extract_ctzt_lower(list(sample_data.loc[expression_table.columns]['time point'])),
        "seq": "PolyA",
        "sex": "F",
        "light": "LD",
        "age_low": 21,
        "age_high": 21,
    },

    "Mezhnina22": {
        "GSE": "GSE211975",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.diet == 'AL'),
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "seq": "PolyA",
        "sex": "F",
        "light": "LD",
        "age_low": 21,
        "age_high": 21,
    },

    "Koritala22": {
        "GSE": "GSE214530",
        "tissue": "Liver",
        "sample_selector": lambda x: (x.treatment == 'Normoxia') and (x.tissue == 'Liver'),
        #"time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])), unknown!!
        # TODO: author contacted, needs to be updated
    },

    "Koronowski22": {
        "GSE": "GSE158600",
        "tissue": "Liver",
        "sample_selector": lambda x: (x['genotype/variation'] == 'WT'),
        "time": lambda sample_data, expression_table: extract_TP(list(sample_data.loc[expression_table.columns]['time of harvest'])),
        "seq": "RiboZero",
        "sex": "M",
        "light": "LD",
        "age_low": 8,
        "age_high": 12,
    },

    "Petrus22": {
        "GSE": "GSE196430",
        "tissue": "Liver",
        "sample_selector": lambda x: (x['genotype/variation'] == 'WT'),
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "seq": "PolyA",
        "sex":  "M",
        "age_low": 8,
        "age_high": 12,
        #TODO: DD? Author contacted.
    },

    "Mekbib22A": {
        "GSE": "GSE182834",
        "tissue": "Liver",
        "sample_selector": lambda x: (x['genotype'] == 'wild type') and (x.Sex == 'Male'),
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "seq": "PolyA",
        "sex": "M",
        "age_low": 20,
        "age_high": 20,
    },

    "Mekbib22A": {
        "GSE": "GSE182834",
        "tissue": "Liver",
        "sample_selector": lambda x: (x['genotype'] == 'wild type') and (x.Sex == 'Female'),
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time point'])),
        "seq": "PolyA",
        "sex": "F",
        "age_low": 20,
        "age_high": 20,
        "light": "LD",
    },

    "Abe22": {
        "GSE": "GSE199061",
        "tissue": "Liver",
        "sample_selector": lambda x: (x['genotype'] == 'WT'),
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "seq": "PolyA",
        "sex": "M",
        "age_low": 8,
        "age_high": 24,
        "light": "DD",
    },

    "Jouffe22": {
        "GSE": "GSE190221",
        "tissue": "Liver",
        "seq": "RiboZero",
        # TODO: can we separate out the sexes?? no....
    },



    #### Kidney Studies ####

    "Yeung17": {
        "GSE": "GSE100457",
        "GSE_postfix": "-GPL17021",
        "sample_selector": lambda x: x.genotype == "WT" and "RNASeq" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['time'])),
        "tissue": "Kidney",
        "seq": "unknown",
    },

    "Zhang14_RNAseq_Kidney_M": {
        "GSE": "GSE54651",
        "sample_selector": lambda x: x.tissue == "kidney",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Kidney",
        "seq":  "unkown",
    },

    "Castelo-Szekely17": {
        "GSE": "GSE81283",
        "sample_selector": lambda x: "total RNA" in x.title,
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns].title)),
        "tissue": "Kidney",
        "seq":  "unkown",
        "trim_adaptor": "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --match-read-wildcards -m 6",
    },

    "Mermet18_Kidney_NightFeed": {
        "GSE": "GSE101423",
        "GSE_postfix": "-GPL17021",
        "sample_selector":  lambda x: x.genotype == "WT" and x.tissue == "Kidney" and x.source_name_ch1 == "PolyA-selected RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Kidney",
        "seq":  "unkown",
    },

    "Mermet18_Kidney": {
        "GSE": "GSE101423",
        "GSE_postfix": "-GPL19057",
        "sample_selector":  lambda x: x.genotype == "Cry1IntronWT" and x.tissue == "Kidney" and x.source_name_ch1 == "PolyA-selected RNA",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['title'])),
        "tissue": "Kidney",
        "seq":  "unkown",
    },

    "Xin21_Kidney_NightFeed": {
        "GSE": "GSE151221",
        "sample_selector": lambda x: x['dietary regimen'] == "NRF",
        "time": lambda sample_data, expression_table: extract_ctzt(list(sample_data.loc[expression_table.columns]['zeitgeber time in hours'])),
        "tissue": "Kidney",
        "seq":  "RiboZero",
    },
}

for name, target in targets.items():
    # Give default values
    if 'short_name' not in target:
        target['short_name'] = name
    if 'sex' not in target:
        target['sex'] = 'unknown'
    if 'light' not in target:
        target['light'] = 'unknown'
    if 'seq' not in target:
        target['seq'] = 'unknown'
    if 'age_low' not in target:
        target['age_low'] = float("NaN")
    if 'age_high' not in target:
        target['age_high'] = float("NaN")
    if 'seq' not in target:
        taret['seq'] = 'unknown'

# List of studies to perform
studies = [
    # Highlighted studies -- all consistent design
    "Weger18_Liver_M", "Morton20_Liver", "Atger15_AdLib", "Yang16B_M", "Weger20_Bmal1WT", "Weger20_HlfDbpTefWT", "Weger20_CryWT", "Fader19", "Greenwell19_AdLib", "Levine20", "Quagliarini19_NormalDiet", "Quagliarini19_NormalDiet_WTvsKO", "Stubblefield18", "Menet12", "Guan20_Liver", "Hidalgo19",
    # Possible highlights
    "Du14", "Pan19",
    #"Janich15", # Different strain, but still C57BL/6; bad sequencing?
    # Remaining
    "Manella21_Liver", "Lahens15", "Weger18_Liver_F", "Zhang14_RNAseq_Liver_M", "Atger15_NightFeed", "Yang16A_M", "Yang16B_F", "Koronowski19_F", "Xin21_Liver_NightFeed", "Yang20", "Kinouchi18_Liver", "Mermet18", "Benegiamo18", "Cajan16", "Chaix19_AdLib_HFD", "Chaix19_NightFeed_HFD", "Chen19", "Gaucher19_Chronic_Cntrl", "Greenwell19_NightFeed", "Hirako18", "Quagliarini19_HFD", "Quagliarini19_HFD_WTvsKO", "Wu19", "Li19_Young", "Li19_Old",
    # new:
    "Aviram21", "Frazier22A", "Frazier22B", "Wu23",
    "Acosta-Rodriguez22A", "Acosta-Rodriguez22B", "Katsioudi22", "Astafev23", "Mezhnina22", "Koronowski22", "Petrus22", "Mekbib22A", "Mekbib22A", "Abe22",
    # Kidney:
    "Yeung17", "Zhang14_RNAseq_Kidney_M", "Castelo-Szekely17", "Mermet18_Kidney_NightFeed", "Mermet18_Kidney", "Xin21_Kidney_NightFeed",
    ]

# List of all tissues being run
tissues = list(set(targets[study]['tissue'] for study in studies))

def studies_by_tissue(tissue):
    return [study for study in studies if targets[study]['tissue'] == tissue]

def select_tissue(studies):
    def f(wildcards, output):
        return [study for study in studies if targets[study]['tissue'] == wildcards.tissue]
    return f
