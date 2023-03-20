import pathlib

import numpy as np
import pandas
import matplotlib
import pylab

JUST_LEFT_COLOR = "#a53131"
JUST_BOTTOM_COLOR = "#0c1860"
BOTH_COLOR = "grey"

from studies import targets as study_info

import studies

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

# Load
summary = pandas.read_csv(snakemake.input.summary, sep="\t")
summary_swapped = summary.copy()
summary_swapped[['study1', 'study2']] = summary_swapped[['study2', 'study1']]
summary = pandas.concat([
    summary,
    summary_swapped,
])


# Reformat
wide = summary.pivot(index=["study1", "study2"], columns="category", values="counts").fillna(0)
both_rhythmic = wide['change'] + wide['same']
one_rhythmic = wide['gain'] + wide['loss']
total = both_rhythmic + one_rhythmic

all_studies = set(summary.study1.unique()).union(summary.study2.unique())
study_order = {study: i for i,study in enumerate(sorted(all_studies, key=lambda x: studies.targets[x]['short_name']))}

# Load JTK and calculate overlaps
jtk = pandas.read_csv(snakemake.input.all_jtk, sep="\t")
jtk['significant'] = jtk['qvalue'] < 0.05
jtk_overlaps = []
for study1, data1 in jtk.groupby("study"):
    sig1 = set(data1.query('significant').ID)
    for study2, data2 in jtk.groupby("study"):
        if study1 == study2:
            continue
        sig2 = set(data2.query('significant').ID)
        jtk_overlaps.append({
            "study1": study1,
            "study2": study2,
            "both": len(sig1.intersection(sig2)),
            "just1": len(sig1.difference(sig2)),
            "just2": len(sig2.difference(sig1)),
        })
jtk_overlaps = pandas.DataFrame(jtk_overlaps).set_index(['study1', 'study2'])
jtk_overlaps['total'] = jtk_overlaps['both'] + jtk_overlaps['just1'] + jtk_overlaps['just2']
print(jtk_overlaps.head())
# Compute the scaling
jtk_scale = max(total.max(), jtk_overlaps['total'].max())

# Load BooteJTK and calculate overlaps
bootejtk = pandas.read_csv(snakemake.input.all_bootejtk, sep="\t")
bootejtk['significant'] = bootejtk['GammaBH'] < 0.05
bootejtk_overlaps = []
for study1, data1 in bootejtk.groupby("study"):
    sig1 = set(data1.query('significant').ID)
    for study2, data2 in bootejtk.groupby("study"):
        if study1 == study2:
            continue
        sig2 = set(data2.query('significant').ID)
        bootejtk_overlaps.append({
            "study1": study1,
            "study2": study2,
            "both": len(sig1.intersection(sig2)),
            "just1": len(sig1.difference(sig2)),
            "just2": len(sig2.difference(sig1)),
        })
bootejtk_overlaps = pandas.DataFrame(bootejtk_overlaps).set_index(['study1', 'study2'])
bootejtk_overlaps['total'] = bootejtk_overlaps['both'] + bootejtk_overlaps['just1'] + bootejtk_overlaps['just2']
print(bootejtk_overlaps.head())
# Compute the scaling
bootejtk_scale = max(total.max(), bootejtk_overlaps['total'].max())

METHOD_OVERLAPS = {
    "JTK": jtk_overlaps,
    "BooteJTK": bootejtk_overlaps,
}

METHOD_SCALE = {
    "JTK": jtk_scale,
    "BooteJTK": bootejtk_scale,
}


def plot_overlaps(ax, x,y, just1, just2, both, scale):
    # Draw a Venn diagram as rectangles
    # with fractions colored to show overlap
    tot = just1 + just2 + both
    square_size = np.sqrt(tot / scale)
    both_size = (both / tot)
    just1_size = (just1 / tot)
    just2_size = (just2 / tot)
    corner = (1-square_size) / 2
    total_rect = matplotlib.patches.Rectangle(
        [i + corner, j + corner],
        width = square_size,
        height = square_size,
        facecolor = BOTH_COLOR,
    )

    if just2_size < 0.2:
        just2_width = min(1 - just1_size, np.sqrt(just2_size))
        just2_height = just2_size / just2_width
    else:
        just2_width = just2_size
        just2_height = 1
    just2_rect = matplotlib.patches.Rectangle(
        [i + corner, j + corner],
        width = just2_width * square_size,
        height = just2_height * square_size,
        facecolor = JUST_LEFT_COLOR,
    )
    #both_rect = matplotlib.patches.Rectangle(
    #    [i + corner + just2_size, j + corner],
    #    width = both_size,
    #    height = square_size,
    #    facecolor = "black",
    #)
    if just1_size < 0.2:
        just1_width = min(1 - just2_size, np.sqrt(just1_size))
        just1_height = (just1_size / just1_width)
    else:
        just1_width = just1_size
        just1_height = 1
    just1_rect = matplotlib.patches.Rectangle(
        [i + corner + square_size - just1_width * square_size, j + corner],
        width = just1_width * square_size,
        height = just1_height * square_size,
        facecolor = JUST_BOTTOM_COLOR,
    )
    ax.add_patch(total_rect)
    #ax.add_patch(both_rect)
    ax.add_patch(just1_rect)
    ax.add_patch(just2_rect)

for method in ['JTK', 'BooteJTK']:
    scale = METHOD_SCALE[method]
    # Plot - showing squares of mini 'venn diagrams'
    LEGEND_GAP = 1
    LEGEND_SPACE = 4
    fig, ax = pylab.subplots(figsize=(9,9), constrained_layout = True)
    for (study1, study2), data in wide.iterrows():
        i = study_order[study1]
        j = study_order[study2]
        if i > j:
            # Plot the compareRhythms overlaps below diagonal
            tot = total.loc[(study1, study2)]
            both = both_rhythmic.loc[(study1, study2)]
            one = one_rhythmic.loc[(study1, study2)]
            just1 = wide['loss'].loc[(study1, study2)]
            just2 = wide['gain'].loc[(study1, study2)]

            plot_overlaps(ax, i, j, just1, just2, both, scale)
        else:
            # Plot the JTK overlaps above diagonal
            df = METHOD_OVERLAPS[method].loc[(study1, study2)]
            both = df['both']
            just1 = df['just1']
            just2 = df['just2']
            plot_overlaps(ax, i, j, just1, just2, both, scale)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_aspect(1.0)
    ax.set_xlim(0, len(study_order) + LEGEND_SPACE)
    ax.set_ylim(0, len(study_order) + 2)
    ax.set_xticks( np.arange(len(study_order))+0.5)
    ax.set_xticklabels(
        labels=[study_info[study]['short_name'] for study in study_order.keys()],
        rotation=90
    )
    ax.set_yticks( np.arange(len(study_order))+0.5)
    ax.set_yticklabels(
        labels=[study_info[study]['short_name'] for study in study_order.keys()],
    )

    # Manually drawn legend
    x = LEGEND_GAP + len(study_order)
    y = (len(study_order) -6*2) / 2
    ax.add_patch(matplotlib.patches.Rectangle([x,y - np.sqrt(500/scale)/2], np.sqrt(500/scale), np.sqrt(500/scale), color="black"))
    ax.add_patch(matplotlib.patches.Rectangle([x,y+2 - np.sqrt(2000/scale)/2], np.sqrt(2000/scale), np.sqrt(2000/scale), color="black"))
    ax.add_patch(matplotlib.patches.Rectangle([x,y+4 - np.sqrt(5000/scale)/2], np.sqrt(5000/scale), np.sqrt(5000/scale), color="black"))
    ax.add_patch(matplotlib.patches.Rectangle([x,y+6 - np.sqrt(5000/scale)/2], np.sqrt(5000/scale), np.sqrt(5000/scale), color=BOTH_COLOR))
    ax.add_patch(matplotlib.patches.Rectangle([x,y+8 - np.sqrt(5000/scale)/2], np.sqrt(5000/scale), np.sqrt(5000/scale), color=JUST_LEFT_COLOR))
    ax.add_patch(matplotlib.patches.Rectangle([x,y+10 - np.sqrt(5000/scale)/2], np.sqrt(5000/scale), np.sqrt(5000/scale), color=JUST_BOTTOM_COLOR))
    ax.annotate("500 genes", [x+0.5,y], xytext = (5,0), textcoords='offset points', va="center")
    ax.annotate("2000 genes", [x+0.5,y+2], xytext = (5,0), textcoords='offset points', va="center")
    ax.annotate("5000 genes", [x+0.5,y+4], xytext = (5,0), textcoords='offset points', va="center")
    ax.annotate("both", [x+0.5,y+6], xytext = (5,0), textcoords='offset points', va="center")
    ax.annotate("just left", [x+0.5,y+8], xytext = (5,0), textcoords='offset points', va="center")
    ax.annotate("just bottom", [x+0.5,y+10], xytext = (5,0), textcoords='offset points', va="center")

    # Label JTK/compareRhythms sections
    ax.annotate(f"{method} rhythmic", [len(study_order), len(study_order)+1], va="center", ha="right")
    ax.annotate("compareRhythms", [len(study_order)+1, len(study_order)], va="top", ha="center", rotation=90)
    fig.savefig(outdir / f"grid.{method}.png", dpi=300)
