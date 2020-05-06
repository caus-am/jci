# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

import os
import json
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os.path
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score



def replace_infs(a):
    # replaces -inf by minimum finite value -1, and +inf by maximum finite value +1
    tmp = np.isfinite(a)
    if len(a[tmp]) == 0:
        amin = 0
        amax = 0
    else:
        amin = min(a[tmp])
        amax = max(a[tmp])
    a[np.isneginf(a)]=amin-1
    a[np.isposinf(a)]=amax+1
    return(a)


algmodesorder=['fisher','asd-obs','asd-pooled','asd-meta','asd-jci0','asd-jci0nt','asd-jci1','asd-jci1nt','asd-jci12','asd-jci12nt','asd-jci123','asd-jci123kt','asd-pikt','asd-jci1-sc','asd-jci123-sc','fci-obs','fci-pooled','fci-meta','fci-jci0','fci-jci1','fci-jci123','fci-jci123r','lcd-mc','lcd-mccon','lcd-mcsct','lcd-mcsctcon','lcd-sc','lcd-sccon','icp-mc','icp-sc','fci-obs-bs','fci-pooled-bs','fci-meta-bs','fci-jci0-bs','fci-jci1-bs','fci-jci123-bs','fci-jci123r-bs','lcd-mc-bs','lcd-mccon-bs','lcd-mcsct-bs','lcd-mcsctcon-bs','lcd-sc-bs','lcd-sccon-bs','icp-mc-bs','icp-sc-bs','cif-mc','cif-mccon','cif-mcsct','cif-mcsctcon','cif-sc','cif-sccon','cif-mc-bs','cif-mccon-bs','cif-mcsct-bs','cif-mcsctcon-bs','cif-sc-bs','cif-sccon-bs']
algmodes_asdjci={'asd-jci123','asd-jci1','asd-jci1nt','asd-jci12','asd-jci12nt','asd-jci0','asd-jci0nt','asd-jci123kt'}
algmodes_asdsc={'asd-jci123-sc','asd-jci1-sc'}
algmodes_asd=algmodes_asdjci | algmodes_asdsc | {'asd-obs','asd-pooled','asd-meta','asd-pikt'}
algmodes_fcijci={'fci-jci123','fci-jci123r','fci-jci1','fci-jci0','fci-jci123-bs','fci-jci123r-bs','fci-jci1-bs','fci-jci0-bs'}
algmodes_fci=algmodes_fcijci | {'fci-obs','fci-pooled','fci-meta','fci-obs-bs','fci-pooled-bs','fci-meta-bs'}
algmodes_lcdmc={'lcd-mc','lcd-mccon','lcd-mcsct','lcd-mcsctcon','lcd-mc-bs','lcd-mccon-bs','lcd-mcsct-bs','lcd-mcsctcon-bs'}
algmodes_lcdsc={'lcd-sc','lcd-sccon','lcd-sc-bs','lcd-sccon-bs'}
algmodes_cif={'cif-mc','cif-mccon','cif-mcsct','cif-mcsctcon','cif-sc','cif-sccon','cif-mc-bs','cif-mccon-bs','cif-mcsct-bs','cif-mcsctcon-bs','cif-sc-bs','cif-sccon-bs'}
algmodes_icp={'icp-sc','icp-mc','icp-sc-bs','icp-mc-bs'}
algmodes_fisher={'fisher'}
algmodes_all=algmodes_asd | algmodes_fci | algmodes_lcdmc | algmodes_lcdsc | algmodes_cif | algmodes_icp | algmodes_fisher
mydpi=150

def translate_algname(alg):
    asdname = 'ASD'
    if alg == 'fisher':
        talg = 'Fisher'
    elif alg == 'asd-obs':
        talg = asdname + '-obs'
    elif alg == 'asd-pooled':
        talg = asdname + '-pooled'
    elif alg == 'asd-meta':
        talg = asdname + '-meta'
    elif alg == 'asd-jci0':
        talg = asdname + '-JCI0'
    elif alg == 'asd-jci0nt':
        talg = asdname + '-JCI0-nt'
    elif alg == 'asd-jci1':
        talg = asdname + '-JCI1'
    elif alg == 'asd-jci1nt':
        talg = asdname + '-JCI1-nt'
    elif alg == 'asd-jci12':
        talg = asdname + '-JCI12'
    elif alg == 'asd-jci12nt':
        talg = asdname + '-JCI12-nt'
    elif alg == 'asd-jci123':
        talg = asdname + '-JCI123'
    elif alg == 'asd-jci123kt':
        talg = asdname + '-JCI123-kt'
    elif alg == 'asd-pikt':
        talg = asdname + '-pikt'
    elif alg == 'asd-jci1-sc':
        talg = asdname + '-JCI1-sc'
    elif alg == 'asd-jci123-sc':
        talg = asdname + '-JCI123-sc'
    elif alg == 'fci-obs':
        talg = 'FCI-obs'
    elif alg == 'fci-pooled':
        talg = 'FCI-pooled'
    elif alg == 'fci-meta':
        talg = 'FCI-meta'
    elif alg == 'fci-jci0':
        talg = 'FCI-JCI0'
    elif alg == 'fci-jci1':
        talg = 'FCI-JCI1'
    elif alg == 'fci-jci123':
        talg = 'FCI-JCI123'
    elif alg == 'fci-jci123r':
        talg = 'FCI-JCI123r'
    elif alg == 'lcd-mc':
        talg = 'LCD-mc'
    elif alg == 'lcd-mc-con':
        talg = 'LCD-mc-con'
    elif alg == 'lcd-mcsct':
        talg = 'LCD-mc-sct'
    elif alg == 'lcd-mcsctcon':
        talg = 'LCD-mc-sct-con'
    elif alg == 'lcd-sc':
        talg = 'LCD-sc'
    elif alg == 'lcd-sccon':
        talg = 'LCD-sc-con'
    elif alg == 'icp-mc':
        talg = 'ICP-mc'
    elif alg == 'icp-sc':
        talg = 'ICP-sc'
    elif alg == 'fci-obs-bs':
        talg = 'FCI-obs-bs'
    elif alg == 'fci-pooled-bs':
        talg = 'FCI-pooled-bs'
    elif alg == 'fci-meta-bs':
        talg = 'FCI-meta-bs'
    elif alg == 'fci-jci0-bs':
        talg = 'FCI-JCI0-bs'
    elif alg == 'fci-jci1-bs':
        talg = 'FCI-JCI1-bs'
    elif alg == 'fci-jci123-bs':
        talg = 'FCI-JCI123-bs'
    elif alg == 'fci-jci123r-bs':
        talg = 'FCI-JCI123r-bs'
    elif alg == 'lcd-mc-bs':
        talg = 'LCD-mc-bs'
    elif alg == 'lcd-mccon-bs':
        talg = 'LCD-mc-con-bs'
    elif alg == 'lcd-mcsct-bs':
        talg = 'LCD-mc-sct-bs'
    elif alg == 'lcd-mcsctcon-bs':
        talg = 'LCD-mc-sct-con-bs'
    elif alg == 'lcd-sc-bs':
        talg = 'LCD-sc-bs'
    elif alg == 'lcd-sccon-bs':
        talg = 'LCD-sc-con-bs'
    elif alg == 'icp-mc-bs':
        talg = 'ICP-mc-bs'
    elif alg == 'icp-sc-bs':
        talg = 'ICP-sc-bs'
    else:
        talg = alg
    return talg

def ordered(algmodes):
    return [s for s in algmodesorder if s in algmodes]

def is_bs(algmode):
    return algmode[-3:] == '-bs'

def has_arel(algmode):
    return algmode in algmodes_all

def has_edge(algmode):
    return algmode in algmodes_asd | algmodes_fci | algmodes_lcdmc | algmodes_lcdsc | algmodes_cif

def has_conf(algmode):
    return algmode in algmodes_asd | algmodes_fci | algmodes_lcdmc | algmodes_lcdsc | algmodes_cif

def has_edgetype(algmode,edgetype):
    if edgetype=='arel':
        return has_arel(algmode)
    elif edgetype=='edge':
        return has_edge(algmode)
    elif edgetype=='conf':
        return has_conf(algmode)
    else:
        return False

def has_rels(algmode,rels):
    if rels=='sys2sys':
        return algmode in algmodes_all - algmodes_fisher
    elif rels=='con2sys':
        return algmode in algmodes_fisher | algmodes_asdjci | algmodes_fcijci | algmodes_lcdmc

def get_order(algmode):
    return algmodesorder.index(algmode)

def get_ls(algmode):
    ls='-'
    if not is_bs(algmode):
        if algmode[0:3] == 'fci':
            ls='-' #'--'
        else:
            ls='-'
    else:
        if algmode[0:3] == 'fci':
            ls='-' #':' #'-.'
        else:
            ls='-' # ':'
    return(ls)

def get_marker(algmode):
    if not is_bs(algmode):
        marker='x'
    else:
        marker='+'
    return(marker)

def get_col(algmode):
    if algmode == 'asd-obs' or algmode == 'fci-obs' or algmode == 'fci-obs-bs':
        col='black'
    elif algmode == 'asd-pooled' or algmode == 'fci-pooled' or algmode == 'fci-pooled-bs':
        col='orange'
    elif algmode == 'asd-meta' or algmode == 'fci-meta' or algmode == 'fci-meta-bs':
        col='magenta'
    elif algmode == 'asd-jci123':
        col='red'
    elif algmode == 'fci-jci123' or algmode == 'fci-jci123-bs':
        col='blue'
    elif algmode == 'fci-jci123r' or algmode == 'fci-jci123r-bs':
        col='green'
    elif algmode == 'asd-jci13':
        col='purple'
    elif algmode == 'asd-jci0' or algmode == 'fci-jci0' or algmode == 'fci-jci0-bs':
        col='lime'
    elif algmode == 'asd-jci0nt':
        col='yellow'
    elif algmode == 'asd-jci1' or algmode == 'fci-jci1' or algmode == 'fci-jci1-bs':
        col='grey'
    elif algmode == 'asd-jci1nt':
        col='cyan'
    elif algmode == 'asd-jci12':
        col='green'
    elif algmode == 'asd-jci12nt':
        col='black'
    elif algmode == 'asd-jci123-sc':
        col='greenyellow'
    elif algmode == 'asd-jci1-sc':
        col='midnightblue'
    elif algmode == 'asd-jci123kt':
        col='brown'
    elif algmode == 'asd-pikt':
        col='blue'
    elif algmode == 'icp-sc' or algmode == 'icp-sc-bs':
        col='green'
    elif algmode == 'icp-mc' or algmode == 'icp-mc-bs':
        col='blue'
    elif algmode == 'fisher':
        col='black'
    elif algmode == 'lcd-mc' or algmode == 'lcd-mc-bs':
        col='black'
    elif algmode == 'lcd-mcsct' or algmode == 'lcd-mcsct-bs':
        col='yellow'
    elif algmode == 'lcd-mccon' or algmode == 'lcd-mccon-bs':
        col='cyan'
    elif algmode == 'lcd-mcsctcon' or algmode == 'lcd-mcsctcon-bs':
        col='grey'
    elif algmode == 'lcd-sc' or algmode == 'lcd-sc-bs':
        col='purple'
    elif algmode == 'lcd-sccon' or algmode == 'lcd-sccon-bs':
        col='red'
    elif algmode == 'cif-mc' or algmode == 'cif-mc-bs':
        col='black'
    elif algmode == 'cif-mcsct' or algmode == 'cif-mcsct-bs':
        col='blue'
    elif algmode == 'cif-mccon' or algmode == 'cif-mccon-bs':
        col='cyan'
    elif algmode == 'cif-mcsctcon' or algmode == 'cif-mcsctcon-bs':
        col='grey'
    elif algmode == 'cif-sc' or algmode == 'cif-sc-bs':
        col='purple'
    elif algmode == 'cif-sccon' or algmode == 'cif-sccon-bs':
        col='red'
    else:
        print('Unknown color %s' % algmode)
    return(col)
    
def get_label(algmode):
    label = algmode
    return(label)


def plot_per_nContext(outdir,plotdir,exppre,exppost,base,only_cyclic,min_pContext,max_pContext,whichOnes,close,texfile,zoom=0,xlim1=[-0.05,1.05],ylim1=[-0.05,1.05],xlim0=[-0.05,1.05],ylim0=[-0.05,1.05]):
    tf = open(texfile,'a')

    for rels in ['sys2sys','con2sys']:
        if rels == 'sys2sys':
            pContexts = range(min_pContext,max_pContext+1)
        elif rels == 'con2sys':
            pContexts = range(max(min_pContext,1),max_pContext+1)
        if rels == 'sys2sys':
            tf.write('\\section{Sys2sys}\n')
        elif rels == 'con2sys':
            tf.write('\\section{Con2sys}\n')
        for pContext in pContexts:
            tf.write('\\subsection{nContext=%d}\n' % pContext)
            exp=exppre+str(pContext)+exppost
            path=outdir+exp+'/'

            if rels == 'sys2sys':
                edge_types = ['arel','edge','conf']
            elif rels == 'con2sys':
                edge_types = ['arel','edge']
            for et in range(len(edge_types)):
                edge_type = edge_types[et]
                if edge_type == 'arel':
                    tf.write('\subsubsection{Ancestral causal relations}\n')
                elif edge_type == 'edge':
                    tf.write('\subsubsection{Direct causes}\n')
                elif edge_type == 'conf':
                    tf.write('\subsubsection{Confounders}\n')
                fignum = 4*len(edge_types)*pContext + 4*et

                for wO in whichOnes:
                    print(pContext,' ',rels,' ',edge_type,' ',wO)
                    algmodes = []
                    if wO == 'asd':
                        algmodes = {s for s in algmodes_asd if has_rels(s,rels)}
                    elif wO == 'asdjci' and rels == 'sys2sys':
                        algmodes = {'asd-obs','asd-pooled','asd-meta','asd-pikt','asd-jci123kt','asd-jci123'}
                    elif wO == 'asdjciN':
                        if rels == 'sys2sys':
                            algmodes = {'asd-obs','asd-jci0','asd-jci1','asd-jci12','asd-jci123','asd-jci123kt'}
                        else:
                            if edge_type == 'edge':
                                algmodes = {'asd-jci0','asd-jci1','asd-jci12','asd-jci123'}
                            else:
                                algmodes = {'fisher','asd-jci0','asd-jci1','asd-jci12','asd-jci123','asd-jci123kt'}
                    elif wO == 'asdjcint':
                        algmodes = {'asd-jci123','asd-jci12','asd-jci12nt','asd-jci1','asd-jci1nt','asd-jci0','asd-jci0nt'}
                    elif wO == 'asdjcisc' and rels == 'sys2sys':
                        algmodes = {'asd-jci1','asd-jci1-sc','asd-jci123','asd-jci123-sc'}
                    elif wO == 'fci':
                        algmodes = {s for s in algmodes_fci if has_rels(s,rels)} #and not is_bs(s)
                        if rels == 'con2sys' and edge_type == 'arel':
                            algmodes.add('fisher')
#                    elif wO == 'fci-bs':
#                        if edge_type == 'arel':
#                            algmodes = {s for s in algmodes_fci if has_rels(s,rels) and is_bs(s)} 
#                        if rels == 'con2sys' and edge_type == 'arel':
#                            algmodes.add('fisher')
                    elif wO == 'lcdicp' and ((rels == 'sys2sys' and edge_type != 'edge') or (rels == 'con2sys' and edge_type != 'arel')):
                        if rels == 'sys2sys':
                            algmodes = {'lcd-sc','lcd-mc','icp-sc','icp-mc'}
                        else:
                            algmodes = {'lcd-mc','icp-mc'}
                    elif wO == 'lcdicp-bs' and ((rels == 'sys2sys' and edge_type != 'edge') or (rels == 'con2sys' and edge_type != 'arel')):
                        if rels == 'sys2sys':
                            algmodes = {'lcd-sc-bs','lcd-mc-bs','icp-sc-bs','icp-mc-bs'}
                        else:
                            algmodes = {'lcd-mc-bs','icp-mc-bs'}
                    elif wO == 'lcd':
                        algmodes = {s for s in algmodes_lcdmc | algmodes_lcdsc | algmodes_icp if has_rels(s,rels) and has_edgetype(s,edge_type) and not is_bs(s)}
                    elif wO == 'lcd-bs':
                        algmodes = {s for s in algmodes_lcdmc | algmodes_lcdsc | algmodes_icp if has_rels(s,rels) and has_edgetype(s,edge_type) and is_bs(s)}
                    elif wO == 'cif':
                        algmodes = {s for s in algmodes_cif | algmodes_icp if has_rels(s,rels) and has_edgetype(s,edge_type) and not is_bs(s)}
                    elif wO == 'cif-bs':
                        algmodes = {s for s in algmodes_cif | algmodes_icp if has_rels(s,rels) and has_edgetype(s,edge_type) and is_bs(s)}
                    elif wO == 'jci':
                        algmodes = {s for s in {'asd-obs','asd-jci123','asd-jci123kt','fci-obs-bs','fci-jci123-bs','fci-jci1-bs','fci-jci0-bs','lcd-mc','icp-sc','lcd-mc-bs','icp-sc-bs','fisher'} if has_rels(s,rels) and has_edgetype(s,edge_type)}
                    elif wO == 'jci1':
                        algmodes = {s for s in {'asd-jci1','fci-jci1-bs','lcd-mc-bs','icp-sc-bs'} if has_rels(s,rels) and has_edgetype(s,edge_type)}
#                    elif wO == 'sc':
#                        algmodes = {s for s in algmodes_sc | algmodes_fisher if has_rels(s,rels) and has_edgetype(s,edge_type) and not is_bs(s)}
#                    elif wO == 'sc-bs':
#                        algmodes = {s for s in algmodes_sc if has_rels(s,rels) and has_edgetype(s,edge_type) and is_bs(s)}
                    elif wO == 'all':
                        algmodes = {s for s in algmodes_all if has_rels(s,rels) and has_edgetype(s,edge_type)}
                    elif wO == 'dit': # direct intervention targets
                        if rels == 'sys2sys':
                            algmodes = {}
                        else:
                            if edge_type == 'edge':
                                algmodes = {'asd-jci0','asd-jci1','asd-jci12','asd-jci123','fci-jci123','fci-jci123-bs'}
                            else:
                                algmodes = {}
                    else:
                        continue

                    anyPlot = 0
                    for algmode in ordered(algmodes):
                        col = get_col(algmode)
                        ls = get_ls(algmode)
                        marker = get_marker(algmode)
                        count,scores,labels,fpr,tpr,rocthr,rocauc,prec0,rec0,pr0thr,avgprec0,prec1,rec1,pr1thr,avgprec1 = roc_pr_stats(path=path,base=base,algmode=algmode,edge_type=edge_type,rels=rels,verbose=0,only_cyclic=only_cyclic)
                        plot_roc_pr(count,scores,labels,fpr,tpr,rocthr,rocauc,prec0,rec0,pr0thr,avgprec0,prec1,rec1,pr1thr,avgprec1,fignum,col,ls,marker=marker,label=('%s' % get_label(algmode)))
                        if count > 0:
                            anyPlot = 1

                    if anyPlot:
                        tf.write(wO+':\\\\\n')
                        tf.write('\\centerline{%\n')
                        decorate_roc(fignum+0,edge_type,rels,plotdir,exp+"_ROC_"+rels+"_"+edge_type+"_"+wO,exp,close,tf)
                        if zoom:
                            if rels == 'sys2sys':
                                decorate_pr(fignum+1,edge_type,rels,1,plotdir,exp+"_PR1_"+rels+"_"+edge_type+"_"+wO+"_zoom",exp,close,tf,xlim1,ylim1)
                                decorate_pr(fignum+2,edge_type,rels,0,plotdir,exp+"_PR0_"+rels+"_"+edge_type+"_"+wO+"_zoom",exp,close,tf,xlim0,ylim0)
                            else:
                                decorate_pr(fignum+1,edge_type,rels,1,plotdir,exp+"_PR1_"+rels+"_"+edge_type+"_"+wO+"_zoom",exp,close,tf,[-0.05,1.05],ylim1)
                                decorate_pr(fignum+2,edge_type,rels,0,plotdir,exp+"_PR0_"+rels+"_"+edge_type+"_"+wO+"_zoom",exp,close,tf,[-0.05,1.05],ylim0)
                        else:
                            decorate_pr(fignum+1,edge_type,rels,1,plotdir,exp+"_PR1_"+rels+"_"+edge_type+"_"+wO,exp,close,tf)
                            decorate_pr(fignum+2,edge_type,rels,0,plotdir,exp+"_PR0_"+rels+"_"+edge_type+"_"+wO,exp,close,tf)
                        draw_legend(fignum+3,plotdir,exp+"_LEG_"+rels+"_"+edge_type+"_"+wO,tf)
                        tf.write('}\n')
                        #plt.figure(fignum)
                        #plt.suptitle('', fontsize=14)
                        #plt.tight_layout()

    tf.close()


def plot_runtimes(outdir,plotdir,exppre,exppost,base,only_cyclic,min_pContext,max_pContext,texfile,algmodes=algmodesorder,label='all',title='Runtimes'):
    fontsize=14
    tf = open(texfile,'a')
    tf.write('\\subsection{Runtimes}\n')
    for pContext in range(min_pContext,max_pContext+1):
        print('Runtimes ',pContext)
        exp=exppre+str(pContext)+exppost
        path=outdir+exp+'/'
        runtimes=get_runtimes(path=path,base=base,only_cyclic=only_cyclic)
        subset = [item for item in runtimes if item in algmodes]
        sorted_subset=sorted(subset,key=get_order,reverse=True)

        fig = plt.figure(pContext,dpi=mydpi)
        #fig.axes[0].set_aspect('equal')
        #plt.xticks([mode for mode in runtimes],rotation='vertical')
        plt.barh(range(len(sorted_subset)),[runtimes[mode] for mode in sorted_subset],tick_label=[translate_algname(mode) for mode in sorted_subset])
        fig.axes[0].set_xscale('log')
        plt.grid(axis='x',which='minor',ls=':', marker='v')
        plt.grid(axis='x',which='major',ls='--', marker='v')
        plt.xlabel('Runtime [s]',fontsize=fontsize)
        plt.ylabel('Method',fontsize=fontsize)
        plt.title(title,fontsize=fontsize)
        plt.yticks(fontsize=7)
        plt.xticks(fontsize=fontsize-2)
        plt.savefig(plotdir+"/pdf/"+exp+"_"+label+"_runtimes.pdf" , transparent=True, format='pdf', bbox_inches="tight", dpi="figure")
        plt.savefig(plotdir+"/png/"+exp+"_"+label+"_runtimes.png" , transparent=True, format='png', bbox_inches="tight", dpi="figure")
        plt.close(pContext)
        tf.write('\\subsubsection{nContext=%d}\n' % pContext)
        tf.write('\\centerline{\\includegraphics[width=\\textwidth]{{'+exp+'_'+label+'_runtimes}.png}}\n')
    tf.close()

def plot_sys2sys_per_algmode(outdir,plotdir,exppre,exppost,base,only_cyclic,min_pContext,max_pContext,close,texfile,xlim=[-0.05,1.05],ylim=[-0.05,1.05]):
    tf = open(texfile,'a')
    tf.write('\\section{As a function of number of context variables}\n')

    for rels in ['sys2sys', 'con2sys']:
        if rels == 'sys2sys':
            edge_types = ['arel','edge','conf']
        elif rels == 'con2sys':
            edge_types = ['arel','edge']
        for et in range(len(edge_types)):
            edge_type = edge_types[et]
            
            algmodes = {s for s in algmodes_all if has_rels(s,rels) and has_edgetype(s,edge_type) and not is_bs(s)}
            if rels == 'sys2sys':
                tf.write('\\subsection{Sys2sys}\n')
            elif rels == 'con2sys':
                tf.write('\\subsection{Con2sys}\n')
            if edge_type == 'arel':
                tf.write('\\subsubsection{Ancestral causal relations}\n') 
            elif edge_type == 'edge':
                tf.write('\\subsubsection{Direct causes}\n') 
            elif edge_type == 'conf':
                tf.write('\\subsubsection{Confounders}\n') 

            for am in ordered(algmodes):
                # fignum = et*len(algmodes)*4 + 4*j
                fignum = 0
                cols = ['purple','red','orange','green','blue','yellow','black','black','black','black','black']
                anyPlot = 0
                print(am,' ',rels,' ',edge_type)
                for k in range(2):
                    if k==0:
                        algmode=am
                    elif k==1:
                        algmode=am+'-bs'
                        
                    if rels == 'sys2sys':
                        pContexts = range(min_pContext,max_pContext+1)
                    elif rels == 'con2sys':
                        pContexts = range(max(min_pContext,1),max_pContext+1)
                    for i in pContexts:
                        exp=exppre+str(i)+exppost
                        path=outdir+exp+'/'
                        col=cols[i]
                        marker = get_marker(algmode)
                        ls = get_ls(algmode)
                        count,scores,labels,fpr,tpr,rocthr,rocauc,prec0,rec0,pr0thr,avgprec0,prec1,rec1,pr1thr,avgprec1 = roc_pr_stats(path=path,base=base,algmode=algmode,edge_type=edge_type,rels=rels,verbose=0,only_cyclic=only_cyclic)
                        plot_roc_pr(count,scores,labels,fpr,tpr,rocthr,rocauc,prec0,rec0,pr0thr,avgprec0,prec1,rec1,pr1thr,avgprec1,fignum,col,ls,marker=marker,label='%d Context Variables' % i)
                        if count > 0:
                            anyPlot = 1

                if anyPlot:
                    tf.write(am+':\\\\\n')
                    tf.write('\\centerline{%\n')
                    decorate_roc(fignum+0,edge_type,rels,plotdir,exppre+"x"+exppost+"_ROC_"+rels+"_"+edge_type+"_"+am,get_label(am),close,tf)
                    decorate_pr(fignum+1,edge_type,rels,1,plotdir,exppre+"x"+exppost+"_PR1_"+rels+"_"+edge_type+"_"+am,get_label(am),close,tf,xlim,ylim)
                    decorate_pr(fignum+2,edge_type,rels,0,plotdir,exppre+"x"+exppost+"_PR0_"+rels+"_"+edge_type+"_"+am,get_label(am),close,tf,xlim,ylim)
                    draw_legend(fignum+3,plotdir,exppre+"x"+exppost+"_LEG_"+rels+"_"+edge_type+"_"+am,tf)
                    tf.write('}\n')
    tf.close()


def plot_runtime_per_algmode(outdir,plotdir,exppre,exppost,base,only_cyclic,min_pContext,max_pContext,texfile):
    tf = open(texfile,'a')
    tf.write('\\subsection{Runtimes}\n')
    fontsize=14
    for algmode in ordered(algmodes_all):
        runtimes={}
        for i in range(min_pContext,max_pContext+1):
            exp=exppre+str(i)+exppost
            path=outdir+exp+'/'
            expruntimes=get_runtimes(path=path,base=base,only_cyclic=only_cyclic,algmode=algmode)
            if algmode in expruntimes:
                runtimes[i]=expruntimes[algmode]
            else:
                runtimes[i]=0
        print(algmode,':',runtimes)

        if len(runtimes) >= 2:
            #plt.xticks([mode for mode in runtimes],rotation='vertical')
            fig = plt.figure(0,dpi=mydpi)
            plt.barh(range(len(runtimes)),[runtimes[i] for i in runtimes],tick_label=[('%d' % i) for i in runtimes])
            plt.grid(axis='x',which='minor',ls=':', marker='v')
            plt.grid(axis='x',which='major',ls='--', marker='v')
            fig.axes[0].set_xscale('log')
            plt.xlabel('Runtime [s]',fontsize=fontsize)
            plt.ylabel('# Context variables',fontsize=fontsize)
            plt.title(translate_algname(algmode),fontsize=fontsize)
            plt.yticks(fontsize=fontsize-2)
            plt.xticks(fontsize=fontsize-2)
            plt.savefig(plotdir+"/pdf/"+exppre+"x"+exppost+"_runtimes_"+algmode+".pdf", transparent=True, format='pdf', bbox_inches="tight", dpi="figure")
            plt.savefig(plotdir+"/png/"+exppre+"x"+exppost+"_runtimes_"+algmode+".png", transparent=True, format='png', bbox_inches="tight", dpi="figure")
            plt.close(0)
            tf.write(algmode+':\\\\\n')
            tf.write('\\centerline{\\includegraphics[width=0.25\\textwidth]{{'+exppre+"x"+exppost+'_runtimes_'+algmode+'}.png}}\n')
    tf.close()


def prcurve(labels,scores,pos_label):
    ind = np.argsort(scores,kind='heapsort')
    cnt = np.count_nonzero(labels==pos_label)
    if cnt == 0:
        rec=np.ones(len(ind))
        prec=np.full(len(ind),fill_value=np.nan)
    else:
        rec=np.zeros(len(ind))
        prec=np.zeros(len(ind))
        thr=np.zeros(len(ind))
        thrcnt = cnt
        j = 0
        for i in range(len(ind)):
            if (i == 0) or (scores[ind[i]] != scores[ind[i-1]]):
                rec[j] = thrcnt * 1.0 / cnt
                prec[j] = thrcnt * 1.0 / (len(ind) - i)
                thr[j] = scores[ind[i]]
                j = j + 1
            if labels[ind[i]] == pos_label:
                thrcnt = thrcnt - 1
        rec=rec[0:j]
        prec=prec[0:j]
        thr=thr[0:j]
    return (prec,rec,thr)

def prcurveold(labels,scores,pos_label):  # should be identical to prcurve, but this one runs in quadratic time
    thr=np.unique(scores) # sorted in ascending order
    rec=np.zeros(len(thr))
    prec=np.zeros(len(thr))
#    rec[0]=len(scores)
#    prec[0]=np.count_nonzero(labels==pos_label)
    for i in range(len(thr)):
        le = len(scores[scores>=thr[i]])
        if le == 0:
            rec[i] = 0
            prec[i] = np.nan
        else:
            cnt = np.count_nonzero(labels==pos_label)
            if cnt == 0:
                rec[i] = 1
                prec[i] = np.nan
            else:
                rec[i] = np.count_nonzero(labels[scores>=thr[i]]==pos_label) / cnt
                prec[i] = np.count_nonzero(labels[scores>=thr[i]]==pos_label) / le
    return (prec,rec,thr)


def roc_pr_stats(path,base,edge_type,algmode,rels,verbose=0,only_cyclic=0):
    count,scores,labels = read_results(path,base,edge_type,algmode,rels,verbose=verbose,only_cyclic=only_cyclic)
#    print('labels:',labels,'scores',scores)
    if count==0:
        return (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    else:
        scores = replace_infs(scores)
#        print('scores',scores)
        fpr,tpr,rocthr = roc_curve(y_true=labels, y_score=scores, pos_label=[1])
        if fpr[0] > 0 or tpr[0] > 0:
            fpr = np.insert(fpr,0,0.0)
            tpr = np.insert(tpr,0,0.0)
            rocthr = np.insert(rocthr,0,np.Inf)
#        fpr = fpr[1:]
#        tpr = tpr[1:]
#        rocthr = rocthr[1:]
        if min(scores)<max(scores):
            rocauc = auc(fpr, tpr)
        else:
            rocauc = np.nan
        #ind1 = scores > 0
        #ind0 = scores < 0
        #prec1, rec1, _ = precision_recall_curve(y_true=labels[ind1], probas_pred=scores[ind1], pos_label=[1])
        #prec0, rec0, _ = precision_recall_curve(y_true=labels[ind0], probas_pred=-scores[ind0], pos_label=[0])
#        Oprec1,Orec1,Opr1thr = precision_recall_curve(y_true=labels, probas_pred=scores, pos_label=[1])
##        pr1thr = np.append(pr1thr,float(max(scores)+1))
#        Oprec1 = Oprec1[:-1]
#        Orec1 = Orec1[:-1]
#        Oprec0,Orec0,Opr0thr = precision_recall_curve(y_true=labels, probas_pred=-scores, pos_label=[0])
##        pr0thr = np.append(pr0thr,float(max(-scores)+1))
#        Oprec0 = Oprec0[:-1]
#        Orec0 = Orec0[:-1]
##        print('ROC:',fpr,tpr,rocthr)
##        print('PR1:',rec1,prec1,pr1thr)
##        print('PR0:',rec0,prec0,pr0thr)
        prec1,rec1,pr1thr=prcurve(labels=labels, scores=scores, pos_label=1)
        prec0,rec0,pr0thr=prcurve(labels=labels, scores=-scores, pos_label=0)
#        if not (np.array_equal(Oprec1,prec1) and np.array_equal(Oprec0,prec0) and np.array_equal(Orec1,rec1) and np.array_equal(Orec0,rec0) and np.array_equal(Opr1thr,pr1thr) and np.array_equal(Opr0thr,pr0thr)):
#            print('Oops!')
#            print('OPR1:',Orec1,Oprec1,Opr1thr)
#            print('OPR0:',Orec0,Oprec0,Opr0thr)
#            print('PR1:',rec1,prec1,pr1thr)
#            print('PR0:',rec0,prec0,pr0thr)
#       NOTE: The difference between the standard and my own implementation is that in case of repeated recall values, the standard method only takes the best precision whereas I just add all
#        prec1O,rec1O,pr1thrO=prcurveold(labels=labels, scores=scores, pos_label=1)
#        prec0O,rec0O,pr0thrO=prcurveold(labels=labels, scores=-scores, pos_label=0)
#        print('PR1:',rec1,prec1,pr1thr)
#        print('PR1old:',rec1O,prec1O,pr1thrO)
#        print('PR1older:',Orec1,Oprec1,Opr1thr)
#        print('PR0:',rec0,prec0,pr0thr)
#        print('PR0old:',rec0O,prec0O,pr0thrO)
#        print('PR0older:',Orec0,Oprec0,Opr0thr)
        avgprec1 = average_precision_score(labels, scores)
        avgprec0 = average_precision_score(labels, -scores)
        return (count,scores,labels,fpr,tpr,rocthr,rocauc,prec0,rec0,pr0thr,avgprec0,prec1,rec1,pr1thr,avgprec1)


def read_results(path,base,edge_type,algmode,rels,verbose=0,only_cyclic=0):
    folders = [folder for folder in os.listdir(path) if re.match('^[0-9]+$',folder)]
    regex = "^"+base+"-"+algmode+"-"+edge_type+".csv$"
    scores = np.empty(shape=0)
    labels = np.empty(shape=0)
    count = 0
    for folder in folders:
        if verbose >= 1:
            print(folder)
        files = [file for file in os.listdir(path+folder) if re.match(regex,file)]
        if len(files) == 1:
            fname = path+folder+'/'+base+".json"
            if verbose >= 2:
                print(fname)
            with open(fname) as json_file:  
                metadata = json.load(json_file)
                systemVars = [x-1 for x in metadata['SystemVars']]
                contextVars = [x-1 for x in metadata['ContextVars']]
                if verbose >= 3:
                    print(systemVars,contextVars)
                
            if only_cyclic != 0:
                fname=path+folder+'/'+base+'-'+'edge'+".csv"
                edges_true=np.loadtxt(fname,delimiter=',',skiprows=1)
                ncc,cclabels=connected_components(edges_true,connection='strong')
                d = edges_true.shape[0]
                if ncc == d: # acyclic
                    if only_cyclic == 1:
                        if verbose >= 1:
                            print('Skipping acyclic model')
                        continue
                    elif only_cyclic == -1:
                        if verbose >= 1:
                            print('Processing acyclic model')
                else:
                    if only_cyclic == 1:
                        if verbose >= 1:
                            print('Processing cyclic model')
                    elif only_cyclic == -1:
                        if verbose >= 1:
                            print('Skipping cyclic model')
                        continue

            fname=path+folder+'/'+base+'-'+edge_type+".csv"
            if verbose >= 2:
                print(fname)
            A_true=np.loadtxt(fname,delimiter=',',skiprows=1)
            if verbose >= 3:
                print(A_true)
            
            fname=path+folder+'/'+base+'-'+algmode+'-'+edge_type+".csv"
            if verbose >= 2:
                print(fname)
            A_pred=np.loadtxt(path+folder+'/'+base+'-'+algmode+'-'+edge_type+".csv",delimiter=',',skiprows=1)
            if verbose >= 3:
                print(A_pred)
            
            if rels=="sys2sys":
                if edge_type == 'conf':
                    N = int(len(systemVars) * (len(systemVars) - 1) / 2)
                else:
                    N = len(systemVars) * (len(systemVars) - 1)
                newlabels = np.full(N,int(0))
                newscores = np.full(N,0.0)
                N = 0
                for row in systemVars:
                    for col in systemVars:
                        if( (edge_type == 'conf' and row < col) or (edge_type != 'conf' and row != col) ):
                            newlabels[N] = int(A_true[row,col])
                            newscores[N] = A_pred[row,col]
                            N+=1
            elif rels=="con2sys":
                N = len(contextVars) * len(systemVars)
                newlabels = np.full(N,int(0))
                newscores = np.full(N,0.0)
                N = 0
                for row in contextVars:
                    for col in systemVars:
                        newlabels[N] = int(A_true[row,col])
                        newscores[N] = A_pred[row,col]
                        N+=1

            labels = np.append(labels,newlabels)
            scores = np.append(scores,newscores)
            count+=1
            if verbose >= 3:
                print(count,labels,scores)
        else:
            if verbose >= 1:
                print('No files found matching %s' % regex)
    return (count,scores,labels)


def get_runtimes(path,base,verbose=0,only_cyclic=0,algmode=''):
    folders = [folder for folder in os.listdir(path) if re.match('^[0-9]+$',folder)]
    if algmode=='':
        regex = "^"+base+"-.*"+".runtime$"
    else:
        regex = "^"+base+"-"+algmode+".runtime$"
    runtimes = {}
    count = 0
    for folder in folders:
        if verbose >= 1:
            print(folder)
        files = [file for file in os.listdir(path+folder) if re.match(regex,file)]
        if len(files) >= 1:
            if only_cyclic != 0:
                fname=path+folder+'/'+base+'-'+'edge'+".csv"
                edges_true=np.loadtxt(fname,delimiter=',',skiprows=1)
                ncc,cclabels=connected_components(edges_true,connection='strong')
                d = edges_true.shape[0]
                if ncc == d: # acyclic
                    if only_cyclic == 1:
                        if verbose >= 1:
                            print('Skipping acyclic model')
                        continue
                    elif only_cyclic == -1:
                        if verbose >= 1:
                            print('Processing acyclic model')
                else:
                    if only_cyclic == 1:
                        if verbose >= 1:
                            print('Processing cyclic model')
                    elif only_cyclic == -1:
                        if verbose >= 1:
                            print('Skipping cyclic model')
                        continue

            for file in files:
                mode=file[len(base)+1:len(file)-8]
                rt=0.0
                with open(path+folder+'/'+file) as f:
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


def plot_roc_pr(count,scores,labels,fpr,tpr,rocthr,rocauc,prec0,rec0,pr0thr,avgprec0,prec1,rec1,pr1thr,avgprec1,fignum,col,ls,marker,label):
    rasterize=False
    markerSize=10
    markerLegend=False
    lineLegend=False
    if count > 0:
        plt.figure(fignum+0,dpi=mydpi)
        # ROC curve
#        print(fpr,tpr,rocthr)
        # start with dotted thin line
        plt.plot(fpr,tpr,color=col,linestyle=':',lw=1,rasterized=rasterize)
        # positives
        if len(rocthr[rocthr>0]) > 0:
            rocthra = min(rocthr[rocthr>0])
            if len(rocthr[rocthr>0]) <= 2:
                plt.plot(fpr[rocthr>=rocthra],tpr[rocthr>=rocthra],color=col,linestyle='None',lw=2,marker=marker,markersize=markerSize,rasterized=rasterize)
                markerLegend = True
            else:
                plt.plot(fpr[rocthr>=rocthra],tpr[rocthr>=rocthra],color=col,linestyle=ls,lw=2,rasterized=rasterize)
                lineLegend = True
        # negatives
        if len(rocthr[rocthr<=0]) > 0:
            rocthrb = max(rocthr[rocthr<=0])
            if len(rocthr[rocthr<=0]) <= 2:
                plt.plot(fpr[rocthr<=rocthrb],tpr[rocthr<=rocthrb],lw=2,color=col,linestyle='None',marker=marker,markersize=markerSize,rasterized=rasterize)
                markerLegend = True
            else:
                plt.plot(fpr[rocthr<=rocthrb],tpr[rocthr<=rocthrb],lw=2,color=col,linestyle=ls,rasterized=rasterize)
                lineLegend = True
        #print (rocthr,min(rocthr[rocthr>0]),max(rocthr[rocthr<0]))
        # PR (positives)
        plt.figure(fignum+1,dpi=mydpi)
#        print(rec1,prec1,pr1thr)
        precbase = (1.0 * len(labels[labels==1])) / len(labels)
        plt.plot([0, 1],[precbase, precbase],color='gray',linestyle=':',label='Random guessing')
        if len(pr1thr[pr1thr>0]) > 0:
            pr1thra = min(pr1thr[pr1thr>0])
            if len(pr1thr[pr1thr>0]) == 1:
                plt.plot(rec1[pr1thr>=pr1thra],prec1[pr1thr>=pr1thra],lw=2,color=col,linestyle='None',alpha=1,marker=marker,markersize=markerSize,rasterized=rasterize)
                markerLegend = True
            else:
                plt.plot(rec1[pr1thr>=pr1thra],prec1[pr1thr>=pr1thra],lw=2,color=col,linestyle=ls,alpha=1,rasterized=rasterize)
                lineLegend = True
        # PR (negatives)
        plt.figure(fignum+2,dpi=mydpi)
#        print(rec0,prec0,pr0thr)
        plt.plot([0, 1],[1-precbase, 1-precbase],color='gray',linestyle=':',label='Random guessing')
        if len(pr0thr[pr0thr>0]) > 0:
            pr0thra = min(pr0thr[pr0thr>0])
            if len(pr0thr[pr0thr>0]) == 1:
                plt.plot(rec0[pr0thr>=pr0thra],prec0[pr0thr>=pr0thra],lw=2,color=col,linestyle='None',alpha=1,marker=marker,markersize=markerSize,rasterized=rasterize)
                markerLegend = True
            else:
                plt.plot(rec0[pr0thr>=pr0thra],prec0[pr0thr>=pr0thra],lw=2,color=col,linestyle=ls,alpha=1,rasterized=rasterize)
                lineLegend = True
        # LEGEND
        plt.figure(fignum+3,dpi=mydpi)
        if markerLegend and lineLegend:
            plt.plot(-1,-1,color=col,linestyle=ls,marker=marker,markersize=markerSize,label=translate_algname(label) + (' (AUC=%0.2f)' % rocauc))
        elif markerLegend and not lineLegend:
            plt.plot(-1,-1,color=col,linestyle='None',marker=marker,markersize=markerSize,label=translate_algname(label) + (' (AUC=%0.2f)' % rocauc))
        elif not markerLegend and lineLegend:
            plt.plot(-1,-1,color=col,linestyle=ls,label=translate_algname(label) + (' (AUC=%0.2f)' % rocauc))

def decorate_roc(num,edge_type,rels,plotdir,basename,title,close,tf):
    fig = plt.figure(num,dpi=mydpi)
    fontsize=14
    if(len(fig.axes) > 0):
        fig.set_size_inches(4,4)
        fig.axes[0].set_aspect('equal')
    #    plt.subplot(143)
        plt.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate',fontsize=fontsize)
        plt.ylabel('True Positive Rate',fontsize=fontsize)
        plt.yticks(fontsize=fontsize-2)
        plt.xticks(fontsize=fontsize-2)
    #    plt.legend(loc="lower right")
        if edge_type == 'arel' and rels == 'sys2sys':
            title='Ancestral causal relations'
        elif edge_type == 'arel' and rels == 'con2sys':
            title='Indirect intervention targets'
        elif edge_type == 'edge' and rels == 'sys2sys':
            title='Direct causal relations'
        elif edge_type == 'edge' and rels == 'con2sys':
            title='Direct intervention targets'
        elif edge_type == 'conf':
            title='Confounders'
        plt.title(title,fontsize=fontsize)
        plt.savefig(plotdir+"/pdf/"+basename+'.pdf', transparent=True, format='pdf', bbox_inches="tight", dpi="figure")
        plt.savefig(plotdir+"/png/"+basename+'.png', transparent=True, format='png', bbox_inches="tight", dpi="figure")
        tf.write('\\includegraphics[width=0.24\\textwidth]{{'+basename+'}.png}\n')
    if close:
        plt.close(num)


def decorate_pr(num,edge_type,rels,pos,plotdir,basename,title,close,tf,xlim=[-0.05,1.05],ylim=[-0.05,1.05]):
    fig = plt.figure(num,dpi=mydpi)
    fontsize=14
    if(len(fig.axes) > 0):
        fig.set_size_inches(4,4)
        if xlim == ylim:
            fig.axes[0].set_aspect('equal')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel('Recall',fontsize=fontsize)
        plt.ylabel('Precision',fontsize=fontsize)
        plt.yticks(fontsize=fontsize-2)
        plt.xticks(fontsize=fontsize-2)
    #    plt.legend(loc="lower right")
    #    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        if edge_type == 'arel' and rels == 'sys2sys':
            if pos:
                title='Ancestors'
            else:
                title='Non-ancestors'
        elif edge_type == 'arel' and rels == 'con2sys':
            if pos:
                title='Indirect intervention targets'
            else:
                title='Indirect intervention non-targets'
        elif edge_type == 'edge' and rels == 'sys2sys':
            if pos:
                title='Direct causes'
            else:
                title='Direct non-causes'
        elif edge_type == 'edge' and rels == 'con2sys':
            if pos:
                title='Direct intervention targets'
            else:
                title='Direct intervention non-targets'
        elif edge_type == 'conf':
            if pos:
                title='Confounders'
            else:
                title='Non-confounders'
        plt.title(title,fontsize=fontsize)
        plt.savefig(plotdir+"/pdf/"+basename+'.pdf', transparent=True, format='pdf', bbox_inches="tight", dpi="figure")
        plt.savefig(plotdir+"/png/"+basename+'.png', transparent=True, format='png', bbox_inches="tight", dpi="figure")
        tf.write('\\includegraphics[width=0.24\\textwidth]{{'+basename+'}.png}\n')
    if close:
        plt.close(num)


def draw_legend(num,plotdir,basename,tf,expand=[-5,-5,5,5]):
    fig = plt.figure(num)
    leg = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
    fig.canvas.draw()
    bbox = leg.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(plotdir+"/pdf/"+basename+'.pdf',transparent=True,format='pdf',dpi="figure",bbox_inches=bbox)
    fig.savefig(plotdir+"/png/"+basename+'.png',transparent=True,format='png',dpi="figure",bbox_inches=bbox)
    tf.write('\\includegraphics[width=0.24\\textwidth]{{'+basename+'}.png}\n')
    plt.close(num)

#def draw_legend(leg,num,edge_type,rels,pos,fname,title):

#   from matplotlib.lines import Line2D
#   custom_lines = [Line2D([0], [0], color='red', lw=4),
#                   Line2D([0], [0], color='blue', lw=4),
#                   Line2D([0], [0], color='green', lw=4)]
#    custom_labels = ['Cold','Medium','Hot']
#
#   fig, ax = plt.figure(num)
#   ax = fig.add_subplot(111)
#
#   leg = ax.legend(custom_lines,custom_labels,frameon=False)
    # hide the axes frame and the x/y labels
#   ax.axis('off')
    #plt.savefig('/tmp/legend.pdf', transparent=True, format='pdf', bbox_inches="tight")
