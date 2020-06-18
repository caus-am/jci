# jci

1. This code belongs to the following paper:

	'Joint Causal Inference from Multiple Contexts'
	by Joris M. Mooij, Sara Magliacane, Tom Claassen
	Journal of Machine Learning Research 21(99):1-108, 2020
	http://jmlr.org/papers/v21/17-123.html

2. Use of this source code is governed by a BSD-style license that can be found
in the LICENSE file. This applies to all code except for that in the file
jci/fci/pcalgjci.R which is licensed as GPL (>= 2) and based upon the R/pcalg.R
file in the source code (pcalg_2.6-10.tar.gz) of the R package pcalg (version 
2.6-10).

3. Eventually we hope to get our extensions of the pcalg package incorporated 
upstream, so that the pcalgjci.R file will become redundant.

4. In order to run this code, you need to 

	a. install the pcalg R package (see more detailed instructions below)

		install.packages('pcalg')

	b. install clingo from https://github.com/potassco/ such that it can 
	be run from the command line

5. Although this contains all code that we used to produce the plots in the
paper, reproducing them is more involved than just running a single command: we
made use of a SLURM cluster and expect that the user will have to make some
changes according to their local setup before everything runs. For example,
one might need to adapt the source code to write the results to another 
directory than the one in /dev/shm/ that is currently hardcoded.

6. This code is released to ensure reproducibility rather than with the aim
of offering user-friendly code to end users. This code might become the basis
of a more user-friendly R package in the future, some day.

7. For questions, please contact Joris M. Mooij <j.m.mooij@uva.nl>


APPENDIX

This is a short example of how to install the R packages in a vanilla
R-base-core 3.5.2 environment in Debian GNU/Linux 10.3. Depending on your OS,
platform, R version, and already installed R packages, you may need to adapt
the procedure. I had to fetch some older versions of some packages (ggm,
glmnet) because the most recent versions only supported R version 3.6 and
higher.

The following sequence of R commands worked for me:

	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install(c('graph','RBGL','cluster','MASS','expm','mgcv','ppcor','Rgraphviz'))

	install.packages('https://cran.r-project.org/src/contrib/Archive/ggm/ggm_2.3.tar.gz', repo=NULL, type='source')
	install.packages('foreach')
	install.packages('https://cran.r-project.org/src/contrib/Archive/glmnet/glmnet_2.0-18.tar.gz', repo=NULL, type='source')

	install.packages(c('pcalg','InvariantCausalPrediction','rjson'))
