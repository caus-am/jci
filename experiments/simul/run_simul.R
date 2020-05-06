# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

args <- commandArgs(trailingOnly = TRUE)

basefilename<-args[1]

pSysObs <- as.numeric(args[2])
pContext <- as.numeric(args[3])
eps <- as.numeric(args[4])
eta <- as.numeric(args[5])
N <- as.numeric(args[6])
acyclic <- as.numeric(args[7])
surgical <- as.numeric(args[8])
seed <- as.numeric(args[9])

cat('run_simul.R: running with arguments',args,'\n')

source('simul_linear_stochastic.R')
simul_linear_stochastic(basefilename,pSysObs,pContext,eps,eta,N,acyclic,surgical,seed)
