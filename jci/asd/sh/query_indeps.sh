#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

if [ "$#" -lt 1 ]; then
	echo 'Usage: query_indeps.sh <p> <inputs...>'
	echo 'Example: query_indeps.sh <p> <inputs...>'
else
	# number of variables
	p=$1
	shift

	TMPFILEOUT=`mktemp /tmp/test.out.XXXXXX`; trap 'rm -f $TMPFILEOUT' 0 15

	pp=`expr $p \- 1`
	for i in `seq 0 $pp`; do
		for j in `seq $[i+1] $pp`; do
			max=$[2**p-1]
			for C in `seq 0 $max`; do
				bla=$[C & (2**i | 2**j)]
				if [ "$bla" -eq "0" ]; then
					M=$[max-2**i-2**j-C]
#					echo running query.sh "indep($i,$j,$C,0,$M,1000000)." $@
#					a=`query.sh "indep($i,$j,$C,0,$M,1000000)." $@` 2> /dev/null
					a=`query.sh ":- th($i,$j,$C,0,$M). :- th($j,$i,$C,0,$M). :- hh($i,$j,$C,0,$M). :- tt($i,$j,$C,0,$M)." $@` 2> /dev/null
#					echo running query.sh "dep($i,$j,$C,0,$M,1000000)." $@
#					b=`query.sh "dep($i,$j,$C,0,$M,1000000)." $@` 2> /dev/null
					b=`query.sh ":- not th($i,$j,$C,0,$M), not th($j,$i,$C,0,$M), not hh($i,$j,$C,0,$M), not tt($i,$j,$C,0,$M)." $@` 2> /dev/null
#					echo $?
					if [ "$a" == "Inf" ]; then
						if [ "$b" == "Inf" ]; then
							echo "% CONTRADICTION: dep($i,$j,$C,0,$M,Inf) AND indep($i,$j,$C,0,$M,Inf)." >> $TMPFILEOUT
						else
							echo -e "% dep($i,$j,$C,0,$M,Inf).\n:- not th($i,$j,$C,0,$M), not th($j,$i,$C,0,$M), not hh($i,$j,$C,0,$M), not tt($i,$j,$C,0,$M)." >> $TMPFILEOUT
						fi
					else
						if [ "$b" == "Inf" ]; then
							echo -e "% indep($i,$j,$C,0,$M,Inf).\n:- th($i,$j,$C,0,$M). :- th($j,$i,$C,0,$M). :- hh($i,$j,$C,0,$M). :- tt($i,$j,$C,0,$M)." >> $TMPFILEOUT
						else
							if [ "$a" -gt "$b" ]; then
								echo "dep($i,$j,$C,0,$M,`expr $a \- $b`)." >> $TMPFILEOUT
							elif [ "$a" -lt "$b" ]; then
								echo "indep($i,$j,$C,0,$M,`expr $b \- $a`)." >> $TMPFILEOUT
							else
								echo "%??dep($i,$j,$C,0,$M)." >> $TMPFILEOUT
							fi
						fi
					fi
				fi
			done
		done
	done

	cat $TMPFILEOUT
	rm -Rf $TMPFILEOUT
fi
