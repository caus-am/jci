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
import os

from plots_commondefs import *

if not os.path.isdir('plots'):
    os.mkdir( 'plots', 0755 );

#for expname in {'ICAM','noICAM','all'}:
for expname in {'noICAM'}:
    make_heatmaps(expname=expname,outpath='out/sachs_' + expname + '/',plotpath='plots/')
