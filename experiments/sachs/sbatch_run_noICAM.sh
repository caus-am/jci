#!/bin/bash
#
# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
#
#SBATCH --job-name=sachs_noICAM
#SBATCH --output=/dev/shm/jmooij1/jci-paper/log/sachs_noICAM-%a.stdout
#SBATCH --error=/dev/shm/jmooij1/jci-paper/log/sachs_noICAM-%a.stderr
#SBATCH --workdir=/dev/shm/jmooij1/jci-paper/code/experiments/sachs
#
#SBATCH --time=30000:00
#SBATCH --mem-per-cpu=1000
#SBATCH --cpus-per-task=1
#SBATCH --exclude=amlab01
#
#SBATCH --mail-type=END,FAIL,REQUEUE,TIME_LIMIT_80
#
#SBATCH --array=0-100

# run as run_sachs_noICAM.sh <miniter> <maxiter>
cd /dev/shm/jmooij1/jci-paper/code/experiments/sachs
./run_sachs_noICAM.sh $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID
