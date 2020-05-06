#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# run as run.sh <pSysObs> <pContext> <eps> <eta> <N> <acyclic> <surgical> <seed> <iters>

echo run.sh called with arguments "$@"

pSysObs=$1
pContext=$2
eps=$3
eta=$4
N=$5
acyclic=$6
surgical=$7
seed=$8
iters=$9

outdir=/dev/shm/jmooij1/jci-paper/out/linear_stochastic_"$pSysObs"_"$pContext"_"$eps"_"$eta"_"$N"_"$acyclic"_"$surgical"/"$seed"
fname='simul'
mkdir -p $outdir

args="$pSysObs $pContext $eps $eta $N $acyclic $surgical"
../Rcmd run_simul.R $outdir/$fname $args $seed

# lcd
for mode in mc mccon mcsct mcsctcon sc sccon; do
  ../Rcmd ../run.R $outdir/$fname lcd $mode 0 $iters
done
# cif
#for mode in mc mccon mcsct mcsctcon sc sccon; do
#  ../Rcmd ../run.R $outdir/$fname cif $mode 0 $iters
#done
# fci
for mode in obs pooled meta jci123 jci1 jci0; do
  ../Rcmd ../run.R $outdir/$fname fci $mode 0 $iters
done
# icp
for mode in sc mc; do
  ../Rcmd ../run.R $outdir/$fname icp $mode 0 $iters
done
# analyze bootstraps
../Rcmd ../analyze.R $outdir/$fname 1
# fisher
../Rcmd ../run.R $outdir/$fname fisher mode
# asd
for mode in obs pooled meta jci123 pikt jci123kt jci13 jci1 jci1nt jci12 jci12nt jci123-sc jci1-sc; do
  ../Rcmd ../run.R $outdir/$fname asd $mode
done
