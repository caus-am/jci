#!/bin/bash
#
# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
#
#SBATCH --job-name=simsl10f
#SBATCH --output=/dev/shm/jmooij1/jci-paper/log/simsl10f-%a.stdout
#SBATCH --error=/dev/shm/jmooij1/jci-paper/log/simsl10f-%a.stderr
#SBATCH --workdir=/dev/shm/jmooij1/jci-paper/code/experiments/simul
#
#SBATCH --time=20000:00
#SBATCH --mem-per-cpu=1000
#SBATCH --cpus-per-task=1
#SBATCH --exclude=amlab01
#
#SBATCH --mail-type=END,FAIL,REQUEUE,TIME_LIMIT_80
#
#SBATCH --array=1-200

# run as run.sh <pSysObs> <pContext> <eps> <eta> <N> <acyclic> <surgical> <seed> <iters>
cd /dev/shm/jmooij1/jci-paper/code/experiments/simul
./run_large.sh 10 10 0.25 0.25 500 1 0 $SLURM_ARRAY_TASK_ID 100
