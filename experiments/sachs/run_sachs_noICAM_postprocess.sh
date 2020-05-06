#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# Use this to postprocess the output of run_sachs_noICAM.sh

fname='sachs'
outdir='out/sachs_noICAM'

# analyze bootstraps
../Rcmd ../analyze.R $outdir/$fname 1
