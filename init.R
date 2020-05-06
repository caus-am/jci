# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# call this by source('init.R',chdir=TRUE)

# JCI preprocessing
source('jci/prepare_jci.R')

# source graphviz functions
source('graphviz/dg2graphviz.R')
source('graphviz/pag2graphviz.R')
source('graphviz/mc2graphviz.R')
source('graphviz/L2graphviz.R')

# source conditional independence tests
source('indepTests/parcortest.R')
source('indepTests/gaussCIFishertest.R')
source('indepTests/gaussCIcontexttest_slow.R')
source('indepTests/gaussCIcontexttest.R')
source('indepTests/gaussCIsincontest.R')

# source jci-fci functions
source('jci/fci_wrapper.R')
source('jci/fci/cpag_to_mc.R')
source('jci/fci/pcalgjci.R')

# source jci-asd functions
source('jci/asd_wrapper.R')
source('jci/asd/asd_indeptests.R')
source('jci/asd/query_clingo.R')
source('jci/asd/run_clingo.R')

# source lcd function
source('jci/lcd.R')

# source icp function
source('jci/icp_wrapper.R')

# source Fisher's method
source('jci/fisher.R')

# define ASP directory
asp_dir<-file.path(getwd(),'jci','asd','asp')
