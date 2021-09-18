import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

import util
from styles import format_study_name

studies = snakemake.params.studies
N_studies = len(studies)

DPI = 150

studies_amp = {}
studies_per = {}
studies_phase = {}
for file, study in zip(snakemake.input.jtk, studies):
    jtk = pandas.read_csv(file, sep="\t", index_col=0)
    p=jtk["ADJ.P"]
    q=jtk["qvalue"]
    amp=jtk["AMP"].copy()
    per=jtk["PER"].copy()
    phase=jtk["LAG"].copy()
    amp[q>0.05]=float("nan")
    per[q>0.05]=float("nan")
    phase[q>0.05]=float("nan")
    studies_amp[study]=amp
    studies_per[study]=per
    studies_phase[study]=phase

def share_xy(axes):
    ''' Force a set of axes to share x and y ranges like pylab.subplots(sharex=True, sharey=True)

    Unlike pylab.subplots, this can be done after-the-fact for performance reasons when you have 
    a large number of axes. Presumably setting sharex=True at the start forces each plot to update
    all existing axes, causing ridiculous extra work.
    '''
    xmin, xmax = float("Inf"), float("-Inf")
    ymin, ymax = float("Inf"), float("-Inf")
    for ax in axes.flatten():
        xmin = min(xmin, ax.get_xlim()[0])
        xmax = max(xmax, ax.get_xlim()[1])
        ymin = min(ymin, ax.get_ylim()[0])
        ymax = max(ymax, ax.get_ylim()[1])
    for ax in axes.flatten():
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

#amplitude
print("Starting amplitude scatter plot")
col=0
fig, axs = pylab.subplots(figsize=(1.5+1.5*N_studies,1.5+1.5*N_studies), nrows=N_studies, ncols=N_studies, squeeze=False)
for study1 in studies:
    amp1 = studies_amp[study1]
    row=0
    for study2 in studies: 
        amp2 = studies_amp[study2]
        axs[col,row].scatter(amp1, amp2, s=1)
        axs[col,row].set_xscale('log')
        axs[col,row].set_yscale('log')
        row+=1
    col+=1
for ax, label in zip(axs[0], studies):
    ax.set_title(format_study_name(label))
for ax, label in zip(axs[:,0], studies): 
    ax.set_ylabel(format_study_name(label))
share_xy(axs)
fig.tight_layout()
fig.savefig(snakemake.output.amplitude, dpi=DPI)

# Not  running now since pretty uninformative: too few distinct values to tell anything
#period
#print("Starting period scatter plot")
#col=0
#fig, axs = pylab.subplots(figsize=(1.5+1.5*N_studies,1.5+1.5*N_studies), nrows=N_studies, ncols=N_studies, squeeze=False)
#for study1 in studies:
#    per1 = studies_per[study1]
#    row=0
#    for study2 in studies:
#        per2 = studies_per[study2]
#        axs[col,row].scatter(per1+numpy.random.normal(size=len(per1))*0.2, per2+numpy.random.normal(size=len(per2))*0.2, s=1)
#        axs[col,row].set_xlim((19, 29))
#        axs[col,row].set_ylim((19, 29))
#        row+=1
#    col+=1
#for ax, label in zip(axs[0], studies):
#    ax.set_title(format_study_name(label))
#for ax, label in zip(axs[:,0], studies):
#    ax.set_ylabel(format_study_name(label))
#share_xy(axs)
#fig.tight_layout()
#fig.savefig(snakemake.output.period, dpi=DPI)

#phase
print("Starting Phase scatter plot")
col=0
fig, axs = pylab.subplots(figsize=(1.5+1.5*N_studies,1.5+1.5*N_studies), nrows=N_studies, ncols=N_studies, squeeze=False)
for study1 in studies:
    phase1 = studies_phase[study1]
    row=0
    for study2 in studies:
        phase2 = studies_phase[study2]
        axs[col,row].scatter(phase1+numpy.random.normal(size=len(phase1))*0.2, phase2+numpy.random.normal(size=len(phase2))*0.2, s=1)
        row+=1
    col+=1
for ax, label in zip(axs[0], studies):
    ax.set_title(format_study_name(label))
for ax, label in zip(axs[:,0], studies):
    ax.set_ylabel(format_study_name(label))
share_xy(axs)
fig.tight_layout()
fig.savefig(snakemake.output.phase, dpi=DPI)
