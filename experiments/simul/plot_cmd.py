# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

import sys
import os
import json
import re
#import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

from plotting import *
rcParams['pdf.fonttype'] = 42 # use Type I fonts for better compatibility

themode = sys.argv[1]

plotdir='plots/'
try:
    os.makedirs(plotdir+'pdf/')
    os.makedirs(plotdir+'png/')
except OSError:
    print('%s exists already' % plotdir)

outdir='out/'
base='simul'
only_cyclic=0

#for themode in ['linear_stochastic_11_9_0.25_0.25_500_1_0','linear_stochastic_11_9_0.15_0.15_500_0_0','linear_stochastic_4_x_0.5_0.5_500_1_1','linear_stochastic_4_x_0.5_0.5_500_1_0','linear_stochastic_4_x_0.5_0.5_500_0_1','linear_stochastic_4_x_0.5_0.5_500_0_0']:
#for themode in ['linear_stochastic_4_x_0.5_0.5_500_1_1','linear_stochastic_4_x_0.5_0.5_500_1_0','linear_stochastic_4_x_0.5_0.5_500_0_1','linear_stochastic_4_x_0.5_0.5_500_0_0']:
#for themode in ['linear_stochastic_11_9_0.25_0.25_500_1_0','linear_stochastic_11_9_0.15_0.15_500_0_0']:
#for themode in ['linear_stochastic_100_10_0.02_0.02_100_1_0']:
#for themode in ['linear_stochastic_4_x_0.5_0.5_500_1_1','linear_stochastic_4_x_0.5_0.5_500_1_0','linear_stochastic_4_x_0.5_0.5_500_0_1','linear_stochastic_4_x_0.5_0.5_500_0_0','linear_stochastic_11_9_0.25_0.25_500_1_0','linear_stochastic_11_9_0.15_0.15_500_0_0','linear_stochastic_11_9_0.25_0.25_500_1_0_gaussCItest','linear_stochastic_11_9_0.15_0.15_500_0_0_gaussCItest','linear_stochastic_100_10_0.02_0.02_100_1_0']:
#for themode in ['linear_stochastic_30_5_0.07_0.07_250_1_0']:
#for themode in ['linear_stochastic_30_5_0.07_0.07_250_1_0_gaussCItest']:
#for themode in ['linear_stochastic_4_x_0.5_0.5_500_1_1']:
#for themode in ['linear_stochastic_4_x_0.5_0.5_500_1_0']:

zoom = 0
if themode == 'linear_stochastic_4_x_0.5_0.5_500_1_1':
    exppre='linear_stochastic_4_'
    exppost='_0.5_0.5_500_1_1'
    min_pContext=0
    max_pContext=4
elif themode == 'linear_stochastic_4_x_0.5_0.5_500_1_0':
    exppre='linear_stochastic_4_'
    exppost='_0.5_0.5_500_1_0'
    min_pContext=0
    max_pContext=4
elif themode == 'linear_stochastic_4_x_0.5_0.5_500_0_1':
    exppre='linear_stochastic_4_'
    exppost='_0.5_0.5_500_0_1'
    min_pContext=0
    max_pContext=4
elif themode == 'linear_stochastic_4_x_0.5_0.5_500_0_0':
    exppre='linear_stochastic_4_'
    exppost='_0.5_0.5_500_0_0'
    min_pContext=0
    max_pContext=4
elif themode == 'linear_stochastic_11_9_0.25_0.25_500_1_0':
    exppre='linear_stochastic_11_'
    exppost='_0.25_0.25_500_1_0'
    min_pContext=9
    max_pContext=9
elif themode == 'linear_stochastic_11_9_0.25_0.25_500_1_0_gaussCItest':
    exppre='linear_stochastic_11_'
    exppost='_0.25_0.25_500_1_0_gaussCItest'
    min_pContext=9
    max_pContext=9
elif themode == 'linear_stochastic_11_9_0.15_0.15_500_0_0':
    exppre='linear_stochastic_11_'
    exppost='_0.15_0.15_500_0_0'
    min_pContext=9
    max_pContext=9
elif themode == 'linear_stochastic_11_9_0.15_0.15_500_0_0_gaussCItest':
    exppre='linear_stochastic_11_'
    exppost='_0.15_0.15_500_0_0_gaussCItest'
    min_pContext=9
    max_pContext=9
elif themode == 'linear_stochastic_10_10_0.25_0.25_500_1_0_gaussCItest':
    exppre='linear_stochastic_10_'
    exppost='_0.25_0.25_500_1_0_gaussCItest'
    min_pContext=10
    max_pContext=10
elif themode == 'linear_stochastic_10_10_0.15_0.15_500_0_0_gaussCItest':
    exppre='linear_stochastic_10_'
    exppost='_0.15_0.15_500_0_0_gaussCItest'
    min_pContext=10
    max_pContext=10
elif themode == 'linear_stochastic_100_10_0.02_0.02_100_1_0':
    exppre='linear_stochastic_100_'
    exppost='_0.02_0.02_100_1_0'
    min_pContext=10
    max_pContext=10
    zoom=1
elif themode == 'linear_stochastic_100_10_0.02_0.02_100_1_0_gaussCItest':
    exppre='linear_stochastic_100_'
    exppost='_0.02_0.02_100_1_0_gaussCItest'
    min_pContext=10
    max_pContext=10
    zoom=1
elif themode == 'linear_stochastic_50_10_0.04_0.04_500_1_0':
    exppre='linear_stochastic_50_'
    exppost='_0.04_0.04_500_1_0'
    min_pContext=10
    max_pContext=10
elif themode == 'linear_stochastic_30_5_0.07_0.07_250_1_0':
    exppre='linear_stochastic_30_'
    exppost='_0.07_0.07_250_1_0'
    min_pContext=5
    max_pContext=5
elif themode == 'linear_stochastic_30_5_0.07_0.07_250_1_0_gaussCItest':
    exppre='linear_stochastic_30_'
    exppost='_0.07_0.07_250_1_0_gaussCItest'
    min_pContext=5
    max_pContext=5
else:
     error('Unknown themode')   
#    elif 0:
#        outdir='out/'
#        exppre='linear_stochastic_4_'
#        exppost='_0.5_0.5_500_1_1_0.0025'
#        base='simul'
#        only_cyclic=0
#        min_pContext=2
#        max_pContext=2
#    elif 0:
#        outdir='out/'
#        exppre='linear_stochastic_3_'
#        exppost='_0.5_0.5_500_0_0'
#        base='simul'
#        only_cyclic=0
#        min_pContext=0
#        max_pContext=3
#    elif 0:
#        outdir='out/'
#        exppre='linear_stochastic_2_'
#        exppost='_0.5_0.5_500_1_1'
#        base='simul'
#        only_cyclic=0
#        min_pContext=0
#        max_pContext=2

print(themode)

texfile=plotdir+exppre+'x'+exppost+'.tex'

tf=open(texfile,'w')
tf.write('\\documentclass{article}\n')
tf.write('\\usepackage{graphicx}\n')
tf.write('\\graphicspath{{png/}}\n')
tf.write('\\begin{document}\n')
tf.close()

#whichOnes = ['dit','fci']
whichOnes = ['asdjci','asdjciN','asdjcisc','asdjcint','lcdicp','lcdicp-bs','fci','dit']
plot_per_nContext(outdir=outdir,plotdir=plotdir,exppre=exppre,exppost=exppost,base=base,only_cyclic=only_cyclic,min_pContext=min_pContext,max_pContext=max_pContext,whichOnes=whichOnes,close=True,texfile=texfile,zoom=0)
if zoom:
    plot_per_nContext(outdir=outdir,plotdir=plotdir,exppre=exppre,exppost=exppost,base=base,only_cyclic=only_cyclic,min_pContext=min_pContext,max_pContext=max_pContext,whichOnes=whichOnes,close=True,texfile=texfile,zoom=zoom,xlim1=[-0.01,0.21],ylim1=[0.74,1.01],xlim0=[-0.05,1.05],ylim0=[0.945,1.005])

algmodes = {'asd-obs','asd-pooled','asd-meta','asd-pikt',
            'asd-jci0','asd-jci1','asd-jci12','asd-jci123','asd-jci123kt','asd-jci1-sc','asd-jci123-sc',
            'fci-obs','fci-pooled','fci-meta','fci-jci0','fci-jci1','fci-jci123',
            'fci-obs-bs','fci-pooled-bs','fci-meta-bs','fci-jci0-bs','fci-jci1-bs','fci-jci123-bs',
            'lcd-mc','lcd-sc','lcd-mc-bs','lcd-sc-bs','icp-mc','icp-sc','icp-mc-bs','icp-sc-bs',
            'fisher'}
plot_runtimes(outdir=outdir,plotdir=plotdir,exppre=exppre,exppost=exppost,base=base,only_cyclic=only_cyclic,min_pContext=min_pContext,max_pContext=max_pContext,texfile=texfile,algmodes=algmodes,label='paper')

plot_sys2sys_per_algmode(outdir=outdir,plotdir=plotdir,exppre=exppre,exppost=exppost,base=base,only_cyclic=only_cyclic,min_pContext=min_pContext,max_pContext=max_pContext,close=True,texfile=texfile)

plot_runtime_per_algmode(outdir=outdir,plotdir=plotdir,exppre=exppre,exppost=exppost,base=base,only_cyclic=only_cyclic,min_pContext=min_pContext,max_pContext=max_pContext,texfile=texfile)

tf=open(texfile,'a')
tf.write('\\end{document}\n')
tf.close()
