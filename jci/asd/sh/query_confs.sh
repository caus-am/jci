#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

if [ "$#" -lt 1 ]; then
	echo 'Usage: query_confs.sh <p> <inputs...>'
	echo 'Example: query_confs.sh <p> <inputs...>'
else
	# number of variables
	p=$1
	shift

	TMPFILEOUT=`mktemp /tmp/test.out.XXXXXX`; trap 'rm -f $TMPFILEOUT' 0 15

	pp=`expr $p \- 1`
	for i in `seq 0 $pp`; do
		for j in `seq $i $pp`; do
			if [ "$i" -ne "$j" ]; then
				a=`query.sh "conf($i,$j)." $@` 2> /dev/null
				b=`query.sh ":-conf($i,$j)." $@` 2> /dev/null
				if [ "$a" == "Inf" ]; then
					if [ "$b" == "Inf" ]; then
						echo "X$i <-> X$j contradiction" >> $TMPFILEOUT
					else
						echo "X$i <-> X$j absent (confidence Inf)" >> $TMPFILEOUT
					fi
				else
					if [ "$b" == "Inf" ]; then
						echo "X$i <-> X$j present (confidence Inf)" >> $TMPFILEOUT
					else
						if [ "$a" -gt "$b" ]; then
							echo "X$i <-> X$j absent (confidence `expr $a \- $b`)" >> $TMPFILEOUT
						elif [ "$a" -lt "$b" ]; then
							echo "X$i <-> X$j present (confidence `expr $b \- $a`)" >> $TMPFILEOUT
						else
							echo "X$i <-> X$j unidentified" >> $TMPFILEOUT
						fi
					fi
				fi
			fi
		done
	done

	cat $TMPFILEOUT
	rm -Rf $TMPFILEOUT
fi
