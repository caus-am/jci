#!/bin/bash
#
# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
#
#SBATCH --job-name=sim4x00
#SBATCH --output=/dev/shm/jmooij1/jci-paper/log/sim4x00-%a.stdout
#SBATCH --error=/dev/shm/jmooij1/jci-paper/log/sim4x00-%a.stderr
#SBATCH --workdir=/dev/shm/jmooij1/jci-paper/code/experiments/simul
#
#SBATCH --time=10000:00
#SBATCH --mem-per-cpu=1000
#SBATCH --cpus-per-task=1
#SBATCH --exclude=amlab01
#
#SBATCH --mail-type=END,FAIL,REQUEUE,TIME_LIMIT_80
#
#SBATCH --array=1-200

# run as run.sh <pSysObs> <pContext> <eps> <eta> <N> <acyclic> <surgical> <seed> <iters>
cd /dev/shm/jmooij1/jci-paper/code/experiments/simul
./run.sh 4 0 0.5 0.5 500 0 0 $SLURM_ARRAY_TASK_ID 100
./run.sh 4 1 0.5 0.5 500 0 0 $SLURM_ARRAY_TASK_ID 100
./run.sh 4 2 0.5 0.5 500 0 0 $SLURM_ARRAY_TASK_ID 100
./run.sh 4 3 0.5 0.5 500 0 0 $SLURM_ARRAY_TASK_ID 100
./run.sh 4 4 0.5 0.5 500 0 0 $SLURM_ARRAY_TASK_ID 100
