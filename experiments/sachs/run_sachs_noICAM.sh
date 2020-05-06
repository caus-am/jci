#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# run as run_sachs_noICAM.sh <miniter> <maxiter>

echo run_sachs_noICAM.sh called with arguments "$@"

datafile=data/sachs-jci-noICAM.csv

fname='sachs'
outdir=/dev/shm/jmooij1/jci-paper/out/sachs_noICAM
mkdir -p $outdir

cp $datafile $outdir/$fname-data.csv
echo -n '{"SystemVars":[1,2,3,4,5,6,7,8,9,10,11],"ContextVars":[13,14,15,16,17,18],"basefilename":"' > $outdir/$fname.json
echo -n "$outdir"/"$fname" >> $outdir/$fname.json
echo '","obsContext":[0,0,0,0,0,0]}' >> $outdir/$fname.json

miniter=$1
maxiter=$2
alpha=1e-2

# lcd
for mode in mc mccon mcsct mcsctcon sc sccon; do
  ../Rcmd ../run.R $outdir/$fname lcd $mode $miniter $maxiter $alpha
done
# cif
for mode in mc mccon mcsct mcsctcon sc sccon; do
  ../Rcmd ../run.R $outdir/$fname cif $mode $miniter $maxiter $alpha
done
## icp
for mode in sc mc; do
  ../Rcmd ../run.R $outdir/$fname icp $mode $miniter $maxiter $alpha
done
# fci
for mode in obs pooled meta jci123 jci1 jci0; do
  ../Rcmd ../run.R $outdir/$fname fci $mode $miniter $maxiter $alpha
done
# fisher
if [ "$1" == "0" ]; then
  ../Rcmd ../run.R $outdir/$fname fisher mode 0 0 $alpha
fi

# analyze bootstraps
../Rcmd ../analyze.R $outdir/$fname 0
# asd
#for mode in obs pooled meta jci123 pikt jci123kt jci13 jci1 jci1nt jci12 jci12nt jci123-sc jci1-sc; do
#  ../Rcmd ../run.R $outdir/$fname asd $mode 0 0 $alpha
#done
