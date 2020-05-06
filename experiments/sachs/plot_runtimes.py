# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

import os
import json
import re
#import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

import sys
sys.path.append('../simul')
from plotting import *
rcParams['pdf.fonttype'] = 42 # use Type I fonts for better compatibility

def plot_sachs_runtimes(outdir,plotdir,exppre,base,algmodes=algmodesorder,label='all',title='Runtimes'):
    exp=exppre
    path=outdir+exp+'/'
    runtimes=get_sachs_runtimes(path=path,base=base)
    print(runtimes)
    sorted_runtimes=sorted(runtimes,key=get_order,reverse=True)

    fig = plt.figure(0,dpi=mydpi)
    #fig.axes[0].set_aspect('equal')
    subset = [item for item in sorted_runtimes if item in algmodes]
    #plt.xticks([mode for mode in runtimes],rotation='vertical')
    plt.barh(range(len(subset)),[runtimes[mode] for mode in subset],tick_label=[translate_algname(mode) for mode in subset])
    fig.axes[0].set_xscale('log')
    plt.grid(axis='x',which='minor',ls=':', marker='v')
    plt.grid(axis='x',which='major',ls='--', marker='v')
    plt.xlabel('Runtime [s]')
    plt.ylabel('Method')
    plt.title(title)
    plt.yticks(fontsize=8)
    plt.savefig(plotdir+exp+"_"+label+"_runtimes.pdf" , transparent=True, format='pdf', bbox_inches="tight", dpi="figure")
    plt.savefig(plotdir+exp+"_"+label+"_runtimes.png" , transparent=True, format='png', bbox_inches="tight", dpi="figure")
    plt.close(0)
    return

def get_sachs_runtimes(path,base,verbose=0,algmode=''):
    if algmode=='':
        regex = "^"+base+"-.*"+".runtime$"
    else:
        regex = "^"+base+"-"+algmode+".runtime$"
    runtimes = {}
    count = 0
    files = [file for file in os.listdir(path) if re.match(regex,file)]
    if len(files) >= 1:
        for file in files:
            mode=file[len(base)+1:len(file)-8]
            rt=0.0
            with open(path+file) as f:
                rt=float(next(f))
                if verbose >= 2:
                    print(rt)
            if( runtimes.get(mode)==None ):
                runtimes[mode]=rt
            else:
                runtimes[mode]+=rt                                

        count+=1
        #print(labels)
        #print(scores)
    for mode in runtimes.keys():
        runtimes[mode]/=count
    return runtimes


plotdir='plots/'
outdir='out/'
base='sachs'

if not os.path.isdir(plotdir):
    os.mkdir(plotdir,0755);

algmodes = {'fci-obs','fci-pooled','fci-meta','fci-jci0','fci-jci1','fci-jci123',
            'fci-obs-bs','fci-pooled-bs','fci-meta-bs','fci-jci0-bs','fci-jci1-bs','fci-jci123-bs',
            'lcd-mc','lcd-sc','lcd-mc-bs','lcd-sc-bs','icp-mc','icp-sc','icp-mc-bs','icp-sc-bs',
            'fisher'}

#for exppre in {'sachs_noICAM','sachs_ICAM','sachs_all'}:
for exppre in {'sachs_noICAM'}:
    plot_sachs_runtimes(outdir=outdir,plotdir=plotdir,exppre=exppre,base=base,algmodes=algmodes,label='paper')
