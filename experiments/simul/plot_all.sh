#!/bin/bash

# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# Use this to produce the plots in the paper

mkdir -p plots

#linear_stochastic_11_9_0.25_0.25_500_1_0 linear_stochastic_11_9_0.15_0.15_500_0_0 linear_stochastic_11_9_0.25_0.25_500_1_0_gaussCItest linear_stochastic_11_9_0.15_0.15_500_0_0_gaussCItest linear_stochastic_100_10_0.02_0.02_100_1_0
for themode in linear_stochastic_4_x_0.5_0.5_500_1_1 linear_stochastic_4_x_0.5_0.5_500_1_0 linear_stochastic_4_x_0.5_0.5_500_0_1 linear_stochastic_4_x_0.5_0.5_500_0_0 linear_stochastic_10_10_0.25_0.25_500_1_0_gaussCItest linear_stochastic_10_10_0.15_0.15_500_0_0_gaussCItest linear_stochastic_100_10_0.02_0.02_100_1_0_gaussCItest; do
	echo $themode
	ipython plot_cmd.py $themode
done
