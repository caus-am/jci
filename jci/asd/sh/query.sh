#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

if [ "$#" -lt 3 ]; then
	echo 'Usage: query.sh <querystr> <input1> <input2> <input3> ...'
	echo 'Example: query.sh "edge(0,1)." 4cycle.ind ../asd_sigma_cyclic.pl'
else
	TMPFILE=`mktemp /tmp/query.XXXXXX`; trap 'rm -f $TMPFILE' 0 15
	TMPFILE2=`mktemp /tmp/query.XXXXXX`; trap 'rm -f $TMPFILE' 0 15
	query=$1
	echo "$query" > $TMPFILE
	shift
	clingo --quiet=2,1 -W no-atom-undefined $TMPFILE $@ > $TMPFILE2
	x=`grep -c UNSATISFIABLE $TMPFILE2`
	y=`grep -c SATISFIABLE $TMPFILE2`
	if [ $x -eq 1 ]; then
		echo Inf
	elif [ $y -eq 1 ]; then
		echo 0
	else
		cat $TMPFILE2 | grep 'Optimization : ' | awk '{print $3}'
	fi
	rm -Rf $TMPFILE2
	rm -Rf $TMPFILE
fi
