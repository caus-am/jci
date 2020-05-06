# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

asd_indeptests <- function(indepfile,data,contextVars,alpha=1e-2,verbose=0,indepTest,suffStat,weightmul=1000,alltests=FALSE,iotests=FALSE) {
  # Prepares for running ASD: does independence tests and writes results to ASP file
  #
  # indepfile:     filename of ASP file with independence test results
  # data:          Nxp matrix
  # contextVars:   indices of context variables (in 1..p)
  # alpha:         p-value threshold for independence test (default: 1e-2)
  # verbose:       verbosity level (default: 0)
  # indepTest:     indepTest function to use
  # suffStat:      sufficient statistics for indepTest
  # weightmul:     weight multiplier for ASP (default: 1000)
  # alltests:      if true, also test (conditional) independences of context variables (default: FALSE)
  # iotests:       if true, only tests (conditional) independences of the form A _||_ B | C with the context variables completely included in the union of B and C

  suppressMessages(library(mgcv))
  
  N<-dim(data)[1]
  stopifnot( N > 0 )
  p<-dim(data)[2]
  stopifnot( p > 0 )

  # run CI tests and write results in ASP format
  if( indepfile != '' ) {
    f<-file(indepfile,'w')
    if( verbose ) {
      cat('Writing independence test results to ',indepfile,'\n',sep='')
    }
  } else
    f<-''
  cat(file=f,'#const nrnodes = ',p,'.\n',sep='')
  if( length(contextVars) > 0 ) {
    for( i in contextVars )
      cat(file=f,'inode(',i-1,').\n',sep='')
  }
  cat(file=f,'% labels: ',colnames(data),'\n')
  cat(file=f,'% ')
  cat(file=f,'\n')
  if( p >= 2 )
    for( i in 1:(p-1) ) {
      for( j in (i+1):p ) {
        if( alltests || !(i %in% contextVars) || !(j %in% contextVars) ) {
          for( Z in 0:(2^p-1) ) {
            if( !bitwAnd(Z,2**(i-1)) && !bitwAnd(Z,2**(j-1)) ) {
              k <- which(intToBits(Z)!=0)
              if( !iotests || (length(setdiff(contextVars,union(c(i,j),k))) == 0) ) { # if iotests, check if all contextvars used in CI statement
		pval<-indepTest(i,j,k,suffStat=suffStat)
		if( !is.na(pval) ) {
		  score <- log(pval) - log(alpha)
		  weight <- round(weightmul * abs(score))
		  if( is.infinite(weight) )
		    weight <- 1000000
		  if( score < 0 ) {
		    cat(file=f,paste('dep(',i-1,',',j-1,',',Z,',0,',2^p-1-2^(i-1)-2^(j-1)-Z,',',as.integer(weight),').\n'))
		  } else {
		    cat(file=f,paste('indep(',i-1,',',j-1,',',Z,',0,',2^p-1-2^(i-1)-2^(j-1)-Z,',',as.integer(weight),').\n'))
		  }
		} else {
		  cat(file=f,paste('%NAdep(',i-1,',',j-1,',',Z,',0,',2^p-1-2^(i-1)-2^(j-1)-Z,',',0,').\n'))
		  warning('NA p-value...')
		}
              }
            }
          }
        }
      }
    }
  if( indepfile != '' )
    close(f)
}


asd_pikt_indeptests <- function(indepfile,data,contextVars,alpha=1e-2,verbose=0,indepTest,suffStat,weightmul=1000,Itargets) {
  # Prepares for running ASD: does independence tests and writes results to ASP file
  #
  # indepfile:     filename of ASP file with independence test results
  # data:          Nxp matrix
  # contextVars:   indices of context variables (in 1..p)
  # alpha:         p-value threshold for independence test (default: 1e-2)
  # verbose:       verbosity level (default: 0)
  # indepTest:     indepTest function to use
  # suffStat:      sufficient statistics for indepTest
  # weightmul:     weight multiplier for ASP (default: 1000)
  # Itargets:      list of intervention targets for each regime

  suppressMessages(library(mgcv))
  
  N<-dim(data)[1]
  stopifnot( N > 0 )
  p<-dim(data)[2]
  stopifnot( p > 0 )

  # run CI tests and write results in ASP format
  if( indepfile != '' ) {
    f<-file(indepfile,'w')
    if( verbose ) {
      cat('Writing independence test results to ',indepfile,'\n',sep='')
    }
  } else
    f<-''
  cat(file=f,'#const nrnodes = ',p,'.\n',sep='')
  if( length(contextVars) > 0 ) {
    for( i in contextVars )
      cat(file=f,'inode(',i-1,').\n',sep='')
  }
  cat(file=f,'% labels: ',colnames(data),'\n')
  cat(file=f,'% ')
  cat(file=f,'\n')
  if( p >= 2 )
    for( i in 1:(p-1) ) {
      for( j in (i+1):p ) {
        if( !(i %in% contextVars) || !(j %in% contextVars) ) {
          for( Z in 0:(2^p-1) ) {
            k <- which(intToBits(Z)!=0)
            if( !(i %in% k) && !(j %in% k) ) {
              nRegimes<-length(suffStat$ns)
              stopifnot(nRegimes >= 1)
              for( R in 1:nRegimes ) {
                pval<-gaussCItest(i,j,k,suffStat=list(C=suffStat$Cs[[R]],n=suffStat$ns[[R]]))
                J<-Itargets[R]
                if( !is.na(pval) ) {
                  score <- log(pval) - log(alpha)
                  weight <- round(weightmul * abs(score))
                  if( is.infinite(weight) )
                    weight <- 1000000
                  if( score < 0 ) {
                    cat(file=f,paste('dep(',i-1,',',j-1,',',Z,',',J,',',2^p-1-2^(i-1)-2^(j-1)-Z,',',as.integer(weight),').\n'))
                  } else {
                    cat(file=f,paste('indep(',i-1,',',j-1,',',Z,',',J,',',2^p-1-2^(i-1)-2^(j-1)-Z,',',as.integer(weight),').\n'))
                  }
                } else {
                  cat(file=f,paste('%NAdep(',i-1,',',j-1,',',Z,',',J,',',2^p-1-2^(i-1)-2^(j-1)-Z,',',0,').\n'))
                  warning('NA p-value...')
                }
              }
            }
          }
        }
      }
    }
  if( indepfile != '' )
    close(f)
}
