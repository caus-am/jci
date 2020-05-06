#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

PATH=$PATH:../jci/asd/sh/

echo "This example corresponds with Fig. 13 in the paper."
echo

echo "JCI-123, acyclic, 3 context variables:"
query_dot.sh 6 chain_3c.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/jci123.pl ../jci/asd/asp/obs_comp_tree.pl
echo
echo "JCI-123, cyclic (sigma-separation), 3 context variables:"
query_dot.sh 6 chain_3c.ind ../jci/asd/asp/asd_sigma_cyclic.pl ../jci/asd/asp/jci123.pl ../jci/asd/asp/obs_comp_tree.pl
echo
echo "JCI-123, acyclic, 1 context variable:"
query_dot.sh 6 chain_1c.ind ../jci/asd/asp/asd_acyclic.pl ../jci/asd/asp/jci123.pl ../jci/asd/asp/obs_comp_tree.pl
echo
echo "JCI-123, cyclic (sigma-separation), 1 context variable:"
query_dot.sh 6 chain_1c.ind ../jci/asd/asp/asd_sigma_cyclic.pl ../jci/asd/asp/jci123.pl ../jci/asd/asp/obs_comp_tree.pl
