import matplotlib
import pylab

from studies import studies, targets

# Determine styles for each study
cmap = pylab.get_cmap("tab20")
marker = ["+", "s", "v", "x"]
linestyles = ["-", "--", "-.", ".."]
color_by_study = {study: cmap(i%20) for i, study in enumerate(studies)}
shape_by_study = {study: marker[i//20] if not targets[study].get('highlight',False) else 'o'
                        for i, study in enumerate(studies)}
linestyle_by_study = {study: linestyles[i//20] for i, study in enumerate(studies)}


# Format study names into multi-lines
def format_study_name(study, max_length=15):
    parts = study.split("_")
    line = parts[0]
    lines = []
    for part in parts[1:]:
        if len(line) + len(part) + 1 > max_length:
            lines.append(line)
            line = part
        else:
            line = line + "_" + part
    lines.append(line)
    joined = '\n'.join(lines)

    if targets[study].get("highlight", False):
        joined = joined + "*"
    return joined
