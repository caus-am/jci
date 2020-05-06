#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

echo "This example corresponds with Fig. 1 in the paper."
echo

PATH=$PATH:../jci/asd/sh/

echo "JCI-1, acyclic:"
query_dot.sh 4 violence.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/jci1.pl ../jci/asd/asp/obs_comp_tree.pl
echo
echo "JCI-1, cyclic (sigma-separation):"
query_dot.sh 4 violence.ind ../jci/asd/asp/asd_sigma_cyclic.pl ../jci/asd/asp/jci1.pl ../jci/asd/asp/obs_comp_tree.pl
