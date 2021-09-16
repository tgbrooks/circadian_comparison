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
