'''
Classifies the type of sequencing performed in a given study
by examining the output from Salmon

I.e. paired or unpaired? Stranded or unstranded?
'''
import json
from studies import studies

metadata = {}
for study, meta_info in zip(studies, snakemake.input.meta_info):
    with open(meta_info) as f:
        metadata[study] = json.load(f)

library_types = {}
study_classification = {}
for study in studies:
    sample_types = [x['library_types'] for x in metadata[study].values()]
    lib_types = list(set(type for type_list in sample_types for type in type_list))
    library_types[study] = lib_types
    if len(lib_types) > 1:
        study_classification[study] = 'mixed'
    else:
        lib_type = lib_types[0]
        classification = ''
        if lib_type[0] in ('I', 'O', 'M'):
            classification = 'paired-'
        else:
            classification = 'unpaired-'
        if 'S' in lib_type:
            classification = classification+'stranded'
        else:
            classification = classification+'unstranded'
        study_classification[study] = classification

with open(snakemake.output[0], "w") as f:
    json.dump(study_classification, f)
