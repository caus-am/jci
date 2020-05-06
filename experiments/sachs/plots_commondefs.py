# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import json
import re
from scipy.linalg import expm
#from matplotlib import rcParams
import os.path
#from scipy.sparse.csgraph import connected_components
#from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

algmodesorder=['fisher','fci-jci123','fci-jci123-bs','fci-jci1','fci-jci1-bs','fci-jci0','fci-jci0-bs','fci-obs','fci-obs-bs','fci-pooled','fci-pooled-bs','fci-meta','fci-meta-bs','lcd-mc','lcd-mc-bs','icp-mc','icp-mc-bs','lcd-sc','lcd-sc-bs','icp-sc','icp-sc-bs','lcd-mccon','lcd-mcsct','lcd-mcsctcon','lcd-sccon','lcd-mccon-bs','lcd-mcsct-bs','lcd-mcsctcon-bs','lcd-sccon-bs','cif-mc','cif-mccon','cif-mcsct','cif-mcsctcon','cif-sc','cif-sccon','cif-mc-bs','cif-mccon-bs','cif-mcsct-bs','cif-mcsctcon-bs','cif-sc-bs','cif-sccon-bs']

def ordered(algmodes):
    return [s for s in algmodesorder if s in algmodes]

def translate_algname(alg):
    if alg == 'em':
        talg = 'Eaton & Murphy (2007)'
    elif alg == 'gt':
        talg = 'Consensus (Sachs et al., 2005)'
    elif alg == 'sachs':
        talg = 'Estimated (Sachs et al., 2005)'
    elif alg == 'mh':
        talg = 'Mooij & Heskes (2013)'
    elif alg == 'mcm':
        talg = 'ACI (Magliacane et al., 2016)'
    elif alg == 'icp_pnas':
        talg = 'ICP (Meinshausen et al., 2016)'
    elif alg == 'hiddenicp_pnas':
        talg = 'hiddenICP (Meinshausen et al., 2016)'
    elif alg == 'fisher':
        talg = 'Fisher'
    elif alg == 'asd-obs':
        talg = 'ASD-obs'
    elif alg == 'asd-pooled':
        talg = 'ASD-pooled'
    elif alg == 'asd-meta':
        talg = 'ASD-meta'
    elif alg == 'asd-jci0':
        talg = 'ASD-JCI0'
    elif alg == 'asd-jci0nt':
        talg = 'ASD-JCI0-nt'
    elif alg == 'asd-jci1':
        talg = 'ASD-JCI1'
    elif alg == 'asd-jci1nt':
        talg = 'ASD-JCI1-nt'
    elif alg == 'asd-jci12':
        talg = 'ASD-JCI12'
    elif alg == 'asd-jci12nt':
        talg = 'ASD-JCI12-nt'
    elif alg == 'asd-jci123':
        talg = 'ASD-JCI123'
    elif alg == 'asd-jci123kt':
        talg = 'ASD-JCI123-kt'
    elif alg == 'asd-pikt':
        talg = 'ASD-pikt'
    elif alg == 'asd-jci1-sc':
        talg = 'ASD-JCI1-sc'
    elif alg == 'asd-jci123-sc':
        talg = 'ASD-JCI123-sc'
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

def flatten(A,rels):
    # flattens square matrix, and removes diagonal if rels == 'sys2sys'
    if rels == 'sys2sys':
        assert(A.shape[0] == A.shape[1])
        N = A.shape[0]
        B = np.zeros(shape=(N*(N-1),1))
        t = 0
        for i in range(N):
            for j in range(N):
                if i != j:
                    B[t] = A[i,j]
                    t = t + 1
    elif rels == 'con2sys':
        B = np.transpose(np.array(np.concatenate(A),ndmin=2))
    return(B)

def read_results(path,base,rels,edge_type,allowedalgs,verbose=0):
    #    folders = [folder for folder in os.listdir(path) if re.match('^[0-9]+$',folder)]
    folders = ['.']
    regex = "^"+base+"-.+-"+edge_type+".csv$"
    scores = np.empty(shape=0)
    labels = []
    for folder in folders:
        if verbose >= 1:
            print(folder)
        fname = path+folder+'/'+base+".json"
        if verbose >= 2:
            print(fname)
        with open(fname) as json_file:  
            metadata = json.load(json_file)
            systemVars = [x-1 for x in metadata['SystemVars']]
            contextVars = [x-1 for x in metadata['ContextVars']]
            systemVars = range(len(systemVars))
            contextVars = range(len(systemVars),len(contextVars) + len(systemVars))
            #if verbose >= 3:
            print(systemVars,contextVars)
        files = [file for file in os.listdir(path+folder) if re.match(regex,file)]
        if verbose >= 1:
            print(files)
        allalgs = [fname.replace(base+"-","").replace("-" + edge_type + ".csv","") for fname in files]
        for alg in ordered(allalgs):
            #        for fname in sorted(files):
            fname = base+'-'+alg+'-'+edge_type+'.csv'
            if verbose >= 2:
                print(fname)
            edges=np.loadtxt(path+folder+'/'+fname,delimiter=',',skiprows=1)
            d = edges.shape[0]
            if verbose >= 3:
                print(edges)
            if rels == 'con2sys' and d == (len(systemVars) + len(contextVars)):
                edges = edges[np.ix_(contextVars,systemVars)]
            elif rels == 'sys2sys' and d >= len(systemVars):
                edges = edges[np.ix_(systemVars,systemVars)]
            else:
                continue
            newlabel = fname.replace(base+"-","").replace("-"+edge_type+".csv","")
            if newlabel in allowedalgs:
                labels.append(newlabel)
                if scores.shape[0] == 0:
                    scores = flatten(edges,rels)
                else:
                    scores = np.hstack([scores, flatten(edges,rels)])
    return (labels,scores)

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

def normalize(scores):
    scoresnorm = scores
    for i in range(scores.shape[1]):
        bla = replace_infs(scoresnorm[:,i])
        scoresnorm[:,i] = bla / (max(-min(bla),max(bla)))
    return scoresnorm

def gt2dot(gt_sys2sys,gt_con2sys,labels):
    # ground truth
    p=gt_sys2sys.shape[0]
    q=gt_con2sys.shape[0]
    print('digraph G {')
    print('subgraph cluster_system {')
    print('label="system";')
    print('style="invis";')
    for i in range(p):
        print(str(i+1) + '[label="' + labels[i] + '"];')
    for i in range(p):
        for j in range(p):
            if gt_sys2sys[i,j]:
                print(str(i+1) + '->' + str(j+1) + ';')
    print('}')
    for k in range(q):
        print(str(k+p+1) + '[label="' + labels[k+p] + '",shape=rectangle];')
    for k in range(q):
        for i in range(p):
            if gt_con2sys[k,i]:
                print(str(k+p+1) + '->' + str(i+1) + ';')
    print('}')
    return

def get_exp(expname):
    if expname=='noICAM':
        exp={'nSys':11,'nCon':6,'labels':["Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK","AKT.inh","G0076","Psitectorigenin","U0126","LY294002","PMA/beta2CAMP + noAlphaCD3/28"]}
    elif expname=='ICAM':
        exp={'nSys':11,'nCon':5,'labels':["Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK","AKT.inh","G0076","Psitectorigenin","U0126","LY294002"]}
    elif expname=='all':
        exp={'nSys':11,'nCon':7,'labels':["Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK","ICAM-2","AKT.inh","G0076","Psitectorigenin","U0126","LY294002","PMA/beta2CAMP + noAlphaCD3/28"]}
    return exp

def get_gt(expname):
    nSys = get_exp(expname)['nSys']
    nCon = get_exp(expname)['nCon']
    labels = get_exp(expname)['labels']
    bl_sys2sys = []
    bl_con2sys = []

    # consensus spp
    bl_sys2sys.append(np.array([[0,1,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,1,0,0],
                                [0,0,0,0,1,0,0,0,1,0,0],
                                [0,0,1,1,0,0,1,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [1,1,0,0,0,1,1,0,0,1,1],
                                [1,1,0,0,0,0,0,0,0,1,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [1,1,0,0,0,0,0,0,0,1,1],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,1,1,0,0,0,0,0,0],
                                [1,1,0,0,0,1,1,0,0,1,1]]))

    # estimated spp
    bl_sys2sys.append(np.array([[0,1,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,1,1,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,1,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [1,1,0,0,0,1,1,0,0,1,1],
                                [1,1,0,0,0,0,0,1,0,1,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))

    # em
    bl_sys2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,1,0,0,0,1,0,0],
                                [0,0,0,0,1,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,1,1,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,1,0,0,0,0,1,0,1,0,1],
                                [0,1,0,0,0,0,0,0,0,1,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,1,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,1,0,1,0,0,0,1],
                                [0,1,1,0,0,1,0,1,0,1,0],
                                [0,0,1,1,0,0,0,0,0,0,0],
                                [1,1,0,0,0,1,1,0,1,0,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [1,0,0,0,1,0,0,1,1,1,0]]))

    # mh
    bl_sys2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [1,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,1,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,1,0,0,0,0,1,0,0,1,1],
                                [1,1,1,1,0,0,1,1,0,1,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [1,1,1,1,0,1,1,1,0,1,1],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [1,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,1,1,0,0,0,0,0,0],
                                [1,1,1,1,0,1,1,1,0,1,1]]))

    # mcm
    bl_sys2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [1,0,0,0,0,1,1,0,0,0,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,1,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,1],
                                [1,1,0,0,0,1,1,0,0,1,1],
                                [1,1,1,1,1,1,1,0,0,1,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,1],
                                [1,1,1,1,1,1,1,0,0,1,1],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [1,0,0,0,0,1,1,0,0,0,1],
                                [0,0,0,1,1,0,0,0,0,0,0],
                                [1,1,1,1,1,1,1,0,0,1,1]]))

    # icp pnas
    bl_sys2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,1,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,1,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,1],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))

    # hiddenicp pnas
    bl_sys2sys.append(np.array([[0,1,0,0,0,0,0,0,0,0,0],
                                [1,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,1,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0],
                                [0,0,0,0,0,0,1,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,1,1],
                                [0,0,0,0,0,0,0,0,1,0,1],
                                [0,0,0,0,0,0,0,0,1,1,0]]))
    bl_con2sys.append(np.array([[0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,0,0,0,0,0,0]]))
        
    # gt em mh mcm icppnas
    bl_complete = []
    bl_complete_anc = []
    bl_sys2sys_anc = []
    bl_con2sys_anc = []
    for T in range(len(bl_sys2sys)):
        if expname=='ICAM':
            bl_con2sys[T] = bl_con2sys[T][1:(nCon+1),0:nSys]
        elif expname=='noICAM':
            bl_con2sys[T] = bl_con2sys[T][1:(nCon+1),0:nSys]
        bl_complete.append(np.concatenate((np.concatenate((bl_sys2sys[T],bl_con2sys[T])),np.zeros(shape=(nSys+nCon,nCon))),axis=1))
        tmp = (expm(bl_complete[T]) > 1e-10) * 1.0
        np.fill_diagonal(tmp,0)
        bl_complete_anc.append(tmp)
        bl_sys2sys_anc.append(tmp[0:nSys,0:nSys])
        bl_con2sys_anc.append(tmp[nSys:(nSys+nCon),0:nSys])
    bl_labels = ['gt','sachs','em','mh','mcm','icp_pnas','hiddenicp_pnas']

    return (bl_sys2sys,bl_sys2sys_anc,bl_con2sys,bl_con2sys_anc,bl_labels)

def get_labels(expname):
    nSys = get_exp(expname)['nSys']
    nCon = get_exp(expname)['nCon']
    labels = get_exp(expname)['labels']

    labels_sys2sys = []
    for i in labels[0:nSys]:
        for j in labels[0:nSys]:
            if i != j:
                labels_sys2sys.append(i + "->" + j)
    labels_con2sys = []
    for i in labels[nSys:(nSys+nCon)]:
        for j in labels[0:nSys]:
            labels_con2sys.append(i + "->" + j)

    return(labels_sys2sys,labels_con2sys)

def make_heatmaps(expname,outpath,plotpath):
    exp=get_exp(expname)
    labels=exp['labels']
    nSys=exp['nSys']
    nCon=exp['nCon']

    (bl_sys2sys,bl_sys2sys_anc,bl_con2sys,bl_con2sys_anc,bl_labels) = get_gt(expname)
    (labels_sys2sys,labels_con2sys) = get_labels(expname)

    modes = [('sys2sys','edge'),('sys2sys','arel'),('con2sys','arel'),('con2sys','edge')]
    for mode in modes:
        rels = mode[0]
        edge_type = mode[1]

        if rels == 'sys2sys':
            if edge_type == 'edge':
                scores_bl = bl_sys2sys
            elif edge_type == 'arel':
                scores_bl = bl_sys2sys_anc
            thelabels = labels_sys2sys
        elif rels == 'con2sys':
            if edge_type == 'edge':
                scores_bl = bl_con2sys
            elif edge_type == 'arel':
                scores_bl = bl_con2sys_anc
            thelabels = labels_con2sys

        if rels == 'sys2sys':
            if edge_type == 'edge':
                allowedalgs = []
                allowedbls = ['gt','sachs','em','mh','icp_pnas','hiddenicp_pnas']
            elif edge_type == 'arel':
                allowedalgs = ['fci-obs','fci-pooled','fci-meta','fci-jci0','fci-jci1','fci-jci123','lcd-mc','lcd-sc','icp-mc','icp-sc','fci-obs-bs','fci-pooled-bs','fci-meta-bs','fci-jci0-bs','fci-jci1-bs','fci-jci123-bs','lcd-mc-bs','lcd-sc-bs','icp-mc-bs','icp-sc-bs']
                allowedbls = ['gt','sachs','em','mh','mcm','icp_pnas','hiddenicp_pnas']
        elif rels == 'con2sys':
            if edge_type == 'edge':
                allowedalgs = ['fci-jci123','fci-jci123-bs']
                allowedbls = ['gt','em']
            elif edge_type == 'arel':
                allowedalgs = ['fisher','fci-jci0','fci-jci1','fci-jci123','fci-jci0-bs','fci-jci1-bs','fci-jci123-bs']
                allowedbls = ['gt','em']
        (algs,scores) = read_results(path=outpath,base='sachs',rels=rels,edge_type=edge_type,allowedalgs=allowedalgs)
        if scores.shape == (0,):
            if rels == 'sys2sys':
                scores = np.empty(shape=(nSys*(nSys-1),0))
            elif rels == 'con2sys':
                scores = np.empty(shape=(nCon*nSys,0))
        for bl in reversed(allowedbls):
            T = bl_labels.index(bl)
            algs.insert(0, bl_labels[T])
            scores = np.insert(scores, 0, values=np.concatenate(flatten(scores_bl[T],rels)), axis=1)

        print(algs)
        print(scores)
        scores = normalize(scores)

        fig, ax = plt.subplots(figsize=(20,20),dpi=200)
        im = ax.imshow(scores,cmap=cm.seismic)
        # We want to show all ticks...
        ax.set_xticks(np.arange(len(algs)))
        ax.set_yticks(np.arange(len(thelabels)))
        # ... and label them with the respective list entries
        if rels == 'sys2sys':
            fs=6
        elif rels == 'con2sys':
            fs=10
        ax.set_xticklabels([translate_algname(alg) for alg in algs],fontsize=fs)
        ax.set_yticklabels(thelabels,fontsize=fs)
        # Rotate the x tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        if rels == 'sys2sys':
            if edge_type == 'arel':
                ax.set_title('Ancestral relations')
            elif edge_type == 'edge':
                ax.set_title('Direct causal relations')
        elif rels == 'con2sys':
            if edge_type == 'arel':
                ax.set_title('Indirect intervention targets')
            elif edge_type == 'edge':
                ax.set_title('Direct intervention targets')
    #    cbar = ax.figure.colorbar(im, ax=ax)
        plt.savefig('plots/sachs_' + expname + '_' + rels + '_' + edge_type + '.pdf', transparent=True, format='pdf', bbox_inches="tight")
        plt.close()

    fig = plt.figure(figsize=(8, 3),dpi=200)
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cm.seismic, norm=mpl.colors.Normalize(vmin=-1, vmax=1), orientation='horizontal')
    cb1.set_label('Confidence Score')
    plt.savefig(plotpath + '/sachs_' + expname + '_LEG.pdf', transparent=True, format='pdf', bbox_inches='tight')
    plt.close()

    gt2dot(bl_sys2sys[0],bl_con2sys[0],labels)

    return
