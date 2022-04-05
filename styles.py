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

color_by_sex = {
        "M": (0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0),
        "F": (1.0, 0.4980392156862745, 0.054901960784313725, 1.0),
        "unknown": (0.5, 0.5, 0.5, 1.0),
}
color_by_light = {
        "LD": (0.8392156862745098, 0.15294117647058825, 0.1568627450980392, 1.0),
        "DD": (0.5803921568627451, 0.403921568627451, 0.7411764705882353, 1.0),
        "unknown": (0.5,0.5,0.5, 1.0),
}

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
