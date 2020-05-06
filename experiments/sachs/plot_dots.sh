#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# Use this to produce the plots in the paper

mkdir -p plots

fname='sachs'
outdir='out/sachs_noICAM'

for y in obs meta pooled; do
	x="$outdir/$fname-fci-$y-pag.dot"
	echo $x
	mv $x $x.bak
	grep -v 'label.*shape' $x.bak > $x
	rm $x.bak
done

for y in jci0 jci1 jci123; do
	sed -i 's/PMA.beta2CAMP...noAlphaCD3.28/PMA\/beta2CAMP + noAlphaCD3\/28/' $outdir/$fname-fci-$y-pag.dot
done

for x in $outdir/$fname*pag.dot; do
	y=`basename $x .dot`
	echo $y
	dot -T eps $x > $y.eps
	epstopdf --outfile=plots/$y.pdf $y.eps
	rm $y.eps
done
