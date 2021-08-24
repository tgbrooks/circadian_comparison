import matplotlib
import pylab

from studies import studies

# Determine styles for each study
cmap = pylab.get_cmap("tab20")
marker = ["o", "+", "s", "v", "x"]
linestyles = ["-", "--", "-.", ".."]
color_by_study = {study: cmap(i%20) for i, study in enumerate(studies)}
shape_by_study = {study: marker[i//20] for i, study in enumerate(studies)}
linestyle_by_study = {study: linestyles[i//20] for i, study in enumerate(studies)}
