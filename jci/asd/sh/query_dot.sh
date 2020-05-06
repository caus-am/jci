#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

if [ "$#" -lt 1 ]; then
	echo 'Usage: query_dot.sh <p> <inputs...>'
	echo 'Example: query_dot.sh <p> <inputs...>'
else
	# number of variables
	p=$1
	shift

	TMPFILEOUT=`mktemp /tmp/test.out.XXXXXX`; trap 'rm -f $TMPFILEOUT' 0 15

	echo 'digraph {' >> $TMPFILEOUT
	grep shape $1 | sed 's/% //' >> $TMPFILEOUT

	pp=`expr $p \- 1`
	for i in `seq 0 $pp`; do
		for j in `seq 0 $pp`; do
			if [ "$i" -ne "$j" ]; then
				a=`query.sh "edge($i,$j)." $@` 2> /dev/null
				b=`query.sh ":-edge($i,$j)." $@` 2> /dev/null
				if [ "$a" == "Inf" ]; then
					if [ "$b" == "Inf" ]; then
#						echo "X$i -> X$j contradiction" >> $TMPFILEOUT
						echo "X$i -> X$j contradiction"
#					else
#						echo "X$i -> X$j absent (confidence Inf)" >> $TMPFILEOUT
					fi
				else
					if [ "$b" == "Inf" ]; then
#						echo "X$i -> X$j present (confidence Inf)" >> $TMPFILEOUT
						echo "$i->$j;" >> $TMPFILEOUT
					else
						if [ "$a" -gt "$b" ]; then
#							echo "X$i -> X$j absent (confidence `expr $a \- $b`)" >> $TMPFILEOUT
#							echo "$i->$j [label=`expr $a \- $b`,color=red];" >> $TMPFILEOUT
							echo -n
						elif [ "$a" -lt "$b" ]; then
#							echo "X$i -> X$j present (confidence `expr $b \- $a`)" >> $TMPFILEOUT
							echo "$i->$j [label=`expr $b \- $a`];" >> $TMPFILEOUT
						else
#							echo "X$i -> X$j unidentified" >> $TMPFILEOUT
							echo "$i->$j [style=\"dashed\"];" >> $TMPFILEOUT
						fi
					fi
				fi
			fi
			if [ "$i" -lt "$j" ]; then
				a=`query.sh "conf($i,$j)." $@` 2> /dev/null
				b=`query.sh ":-conf($i,$j)." $@` 2> /dev/null
				if [ "$a" == "Inf" ]; then
					if [ "$b" == "Inf" ]; then
#						echo "X$i <-> X$j contradiction" >> $TMPFILEOUT
						echo "X$i <-> X$j contradiction"
#					else
#						echo "X$i <-> X$j absent (confidence Inf)" >> $TMPFILEOUT
					fi
				else
					if [ "$b" == "Inf" ]; then
#						echo "X$i <-> X$j present (confidence Inf)" >> $TMPFILEOUT
						echo "$i->$j [dir=\"both\"];" >> $TMPFILEOUT
					else
						if [ "$a" -gt "$b" ]; then
#							echo "X$i <-> X$j absent (confidence `expr $a \- $b`)" >> $TMPFILEOUT
#							echo "$i->$j [dir=\"both\",label=`expr $a \- $b`,color=red];" >> $TMPFILEOUT
							echo -n
						elif [ "$a" -lt "$b" ]; then
#							echo "X$i <-> X$j present (confidence `expr $b \- $a`)" >> $TMPFILEOUT
							echo "$i->$j [dir=\"both\",label=`expr $b \- $a`];" >> $TMPFILEOUT
						else
#							echo "X$i <-> X$j unidentified" >> $TMPFILEOUT
							echo "$i->$j [dir=\"both\",style=\"dashed\"];" >> $TMPFILEOUT
						fi
					fi
				fi
			fi
		done
	done
	
	echo '}' >> $TMPFILEOUT

	cat $TMPFILEOUT
	rm -Rf $TMPFILEOUT
fi
