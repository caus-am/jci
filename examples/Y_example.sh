#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

echo "This example corresponds with Prop. 11 and Fig. 4 in the paper."
echo

PATH=$PATH:../jci/asd/sh/

echo "Acyclic:"
query_dot.sh 4 Y.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/obs_comp_tree.pl
echo
echo "Cyclic (sigma-separation):"
query_dot.sh 4 Y.ind ../jci/asd/asp/asd_sigma_cyclic.pl ../jci/asd/asp/obs_comp_tree.pl
