import json
import pathlib
import pandas
import numpy
import matplotlib
matplotlib.use("Agg")
import pylab

import util

studies = snakemake.params.studies
N_studies = len(studies)

DPI = 150

studies_amp = {}
studies_per = {}
studies_phase = {}
for file, study in zip(snakemake.input.jtk, studies):
    jtk = pandas.read_csv(file, sep="\t", index_col=0)
    p=jtk["ADJ.P"]
    amp=jtk["AMP"]
    per=jtk["PER"]
    phase=jtk["LAG"]
    amp[p>0.05]=float("nan")
    per[p>0.05]=float("nan")
    phase[p>0.05]=float("nan")
    studies_amp[study]=amp
    studies_per[study]=per
    studies_phase[study]=phase

#amplitude
row=0
fig, axs = pylab.subplots(figsize=(1.5+1.5*N_studies,1.5+1.5*N_studies), nrows=N_studies, ncols=N_studies, sharex=True, sharey=True, squeeze=False)
for study1 in studies:
    amp1 = studies_amp[study1]
    col=0
    for study2 in studies: 
        amp2 = studies_amp[study2]
        axs[row,col].scatter(amp1, amp2, s=1)
        axs[row,col].set_xscale('log')
        axs[row,col].set_yscale('log')
        col+=1
    row+=1
for ax, label in zip(axs[0], studies):
    ax.set_title(label)
for ax, label in zip(axs[:,0], studies): 
    ax.set_ylabel(label)
fig.tight_layout()
fig.savefig(snakemake.output.amplitude, dpi=DPI)

#period
row=0
fig, axs = pylab.subplots(figsize=(1.5+1.5*N_studies,1.5+1.5*N_studies), nrows=N_studies, ncols=N_studies, sharex=True, sharey=True, squeeze=False)
for study1 in studies:
    per1 = studies_per[study1]
    col=0
    for study2 in studies:
        per2 = studies_per[study2]
        axs[row,col].scatter(per1+numpy.random.normal(size=len(per1))*0.2, per2+numpy.random.normal(size=len(per2))*0.2, s=1)
        axs[row,col].set_xlim((19, 29))
        axs[row,col].set_ylim((19, 29))
        col+=1
    row+=1
for ax, label in zip(axs[0], studies):
    ax.set_title(label)
for ax, label in zip(axs[:,0], studies):
    ax.set_ylabel(label)
fig.tight_layout()
fig.savefig(snakemake.output.period, dpi=DPI)

#phase
row=0
fig, axs = pylab.subplots(figsize=(1.5+1.5*N_studies,1.5+1.5*N_studies), nrows=N_studies, ncols=N_studies, sharex=True, sharey=True, squeeze=False)
for study1 in studies:
    phase1 = studies_phase[study1]
    col=0
    for study2 in studies:
        phase2 = studies_phase[study2]
        axs[row,col].scatter(phase1+numpy.random.normal(size=len(phase1))*0.2, phase2+numpy.random.normal(size=len(phase2))*0.2, s=1)
        col+=1
    row+=1
for ax, label in zip(axs[0], studies):
    ax.set_title(label)
for ax, label in zip(axs[:,0], studies):
    ax.set_ylabel(label)
fig.tight_layout()
fig.savefig(snakemake.output.phase, dpi=DPI)
