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
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns].time),
    },

    "Weger18": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['zeitgeber time']),
    },

    "Weger18_Liver_M": {
        "GSE": "GSE114400",
        "sample_selector": lambda x: x.Sex == "male" and x['gut microbiome status'] == "conventional raised",
        "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['zeitgeber time']),
    },

    "Weger18_Liver_F": {
    "GSE": "GSE114400",
    "sample_selector": lambda x: x.Sex == "Female" and x['gut microbiome status'] == "conventional raised",
    "time": lambda sample_data, expression_table: list(sample_data.loc[expression_table.columns]['zeitgeber time']),
     },
}

# List of studies to perform
studies = ["Lahens15", "Weger18_Liver_M", "Weger18_Liver_F"]
