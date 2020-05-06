# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

lcd <- function(data,systemVars,contextVars,alpha=1e-2,verbose=0,subsamplefrac=0.0,datamode='merge',test='gaussCIcontexttest',conservative=FALSE) {
  # data:          Nxp matrix
  # systemVars:    indices of sysem variables (in 1..p)
  # contextVars:   indices of context variables (in 1..p)
  # alpha:         p-value threshold for independence test (default: 1e-2)
  # verbose:       verbosity level (default: 0)
  # subsamplefrac: fraction of subsamples used for subsampling; if 0, don't subsample and use all data
  # datamode:      one of {'multiple','merge'}, default='merge'
  # test:          'gaussCItest' / 'gaussCIcontexttest_slow' / 'gaussCIcontexttest' / 'gaussCIsincontest'
  #                   gaussCItest: just partial correlations
  #                   gaussCIcontexttest_slow: tests within each context value separately and combines with Fisher's method
  #                   gaussCIcontexttest: tests within each context value separately and combines with Fisher's method, faster implementation
  #                   gaussCIsincontest: ICP-type test for mixture of Gaussian data
  # conservative:  use conservative p-value? (default: FALSE)

  suppressMessages(library(pcalg))

  # prepare data
  X<-prepare_jci(data,systemVars,contextVars,subsamplefrac,datamode,obsContext=c())
  data<-X$data
  systemVars<-X$systemVars
  contextVars<-X$contextVars
  N<-X$N
  p<-X$p
  removeNAs<-X$removeNAs

  pSys <- length(systemVars)
  pCon <- length(contextVars)
  edge <- matrix(0,p,p)
  colnames(edge) <- colnames(edge)
  arel <- matrix(0,p,p)
  colnames(arel) <- colnames(data)
  conf <- matrix(0,p,p)
  colnames(conf) <- colnames(conf)

  if( pCon > 0 ) {
    # setup independence test
    if( test == 'gaussCItest' ) {
      indepTest<-gaussCItest
      suffStat<-list(C=cor(data),n=N,removeNAs=removeNAs)
    } else if( test == 'gaussCIcontexttest_slow' ) {
      indepTest<-gaussCIcontexttest_slow
      suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
    } else if( test == 'gaussCIcontexttest' ) {
      indepTest<-gaussCIcontexttest
      # find indices for unique joint values of context variables
      uniqueContextValues<-uniquecombs(data[,contextVars])
      regimes<-attr(uniqueContextValues,'index')
      suffStat<-list(data=data,contextVars=contextVars,uniqueContextValues=uniqueContextValues,regimes=regimes,verbose=verbose,removeNAs=removeNAs)
    } else if( test == 'gaussCIsincontest' ) {
      indepTest<-gaussCIsincontest
      suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
    } else {
      error('unknown test')
    }

    # run LCD
    p_ij <- matrix(0,pSys,pSys)
    p_ci <- matrix(0,pCon,pSys)
    for( i in 1:pSys ) {
      for( j in 1:pSys )
        if( i != j )
          p_ij[i,j] <- indepTest(i,j,c(),suffStat)
      for( c in 1:pCon )
        p_ci[c,i] <- indepTest(pSys+c,i,c(),suffStat)
    }
    for( i in 1:pSys )
      for( j in 1:pSys ) 
        if( i != j && p_ij[i,j] < alpha )
          for( c in 1:pCon )
            if( p_ci[c,i] < alpha ) {
              pval <- indepTest(pSys+c,j,c(i),suffStat)
              if( pval >= alpha ) {
#                arel[i,j] <- arel[i,j] + 1
                if( conservative)  {
                  pval2 <- indepTest(pSys+c,i,c(j),suffStat)
                  pval3 <- indepTest(i,j,c(pSys+c),suffStat)
                  newp <- -log(max(c(p_ij[i,j],p_ci[c,i],p_ci[c,j],pval2,pval3)))
                } else {
                  newp <- -log(p_ci[c,j])
                }
                arel[i,j] <- max(arel[i,j],newp)
#               arel[pSys+c,i] <- max(arel[pSys+c,i],newp)  # this only holds under JCI-2
                conf[i,j] <- min(conf[i,j],conf[j,i],-newp)
                conf[j,i] <- conf[i,j]
                edge[pSys+c,j] <- min(edge[pSys+c,j],-newp)
                conf[pSys+c,j] <- min(conf[pSys+c,j],-newp)
                conf[j,pSys+c] <- conf[pSys+c,j]
                if( verbose >= 2 )
                  cat('LCD triple <',pSys+c,',',i,',',j,'> found, pval=',p_ci[c,j],'\n',sep='')
              }
            }
  }

  # return result
  result<-list(p=p,systemVars=systemVars,contextVars=contextVars,labels=colnames(data),arel=arel,conf=conf,edge=edge)
}
