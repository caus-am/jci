#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

PATH=$PATH:../jci/asd/sh/

echo "This example corresponds with Fig. 14 in the paper."
echo

echo "JCI123, acyclic:"
query_dot.sh 3 indc.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/jci123.pl ../jci/asd/asp/obs_comp_tree.pl
echo
echo "JCI12, independent context variables, acyclic:"
query_dot.sh 3 indc.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/jci12indc.pl ../jci/asd/asp/obs_comp_tree.pl

echo
echo "Find all optimals models (JCI12, independent context variables, acyclic):"
clingo --opt-mode=optN indc.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/jci12indc.pl ../jci/asd/asp/obs_comp_tree.pl ../jci/asd/asp/show_minimal.pl
