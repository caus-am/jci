# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

fci_wrapper <- function(data,systemVars,contextVars,alpha=1e-2,verbose=0,subsamplefrac=0.0,mode='obs',test='gaussCIcontexttest',obsContext=matrix(0,1,length(contextVars)),doPDsep=TRUE) {
  # data:          Nxp matrix
  # systemVars:    indices of sysem variables (in 1..p)
  # contextVars:   indices of context variables (in 1..p)
  # alpha:         p-value threshold for independence test (default: 1e-2)
  # verbose:       verbosity level (default: 0)
  # subsamplefrac: fraction of subsamples used for subsampling; if 0, don't subsample and use all data
  # mode:          'obs' / 'pooled' / 'meta' / 'jci123' / 'jci123r' / 'jci1' / 'jci0'
  # test:          'gaussCItest' / 'gaussCIcontexttest_slow' / 'gaussCIcontexttest' (only relevant for mode=='jci*')
  #                   gaussCItest: just partial correlations
  #                   gaussCIcontexttest_slow: tests within each context value separately and combines with Fisher's method
  #                   gaussCIcontexttest: tests within each context value separately and combines with Fisher's method, faster implementation
  # obsContext:    vector of values of context variables that defines the observational context
  # doPDsep:       boolean (whether PDsep phase should be done in FCI; default: TRUE)

  suppressMessages(library(pcalg))
  suppressMessages(library(mgcv))

  # prepare data
  if( mode == 'obs' )
    X<-prepare_jci(data,systemVars,contextVars,subsamplefrac,'obsonly',obsContext)
  else
    X<-prepare_jci(data,systemVars,contextVars,subsamplefrac,'multiple',obsContext=c())
  data<-X$data
  systemVars<-X$systemVars
  contextVars<-X$contextVars
  N<-X$N
  p<-X$p
  removeNAs<-X$removeNAs

  # setup independence test
  if( mode == 'obs' || mode == 'pooled' || length(contextVars) == 0 ) {
    if( removeNAs )
      stop('removeNAs not implemented yet for this mode')
    indepTest<-gaussCItest
    suffStat<-list(C=cor(data[,systemVars]),n=N)
  } else if( mode == 'meta' ) {
    if( removeNAs )
      stop('removeNAs not implemented yet for this mode')
    # find indices for unique joint values of context variables
    uniqueContextValues<-uniquecombs(data[,contextVars])
    regimes<-attr(uniqueContextValues,'index')
    nRegimes<-max(regimes)

    indepTest<-gaussCIFishertest # use Fisher's method for combining p-values tested in different regimes
    Cs<-vector(mode='list',length=nRegimes)
    ns<-vector(mode='integer',length=nRegimes)
    for( R in 1:nRegimes ) {
      subset<-which(regimes==R)
      Cs[[R]]<-cor(data[subset,systemVars])
      ns[[R]]<-length(subset)
    }
    suffStat<-list(Cs=Cs,ns=ns,verbose=verbose)
  } else if( mode == 'jci123' || mode == 'jci123r' || mode == 'jci1' || mode == 'jci0' ) {
    if( mode == 'jci123r' && test != 'gaussCIcontexttest' )
      stop('not implemented yet')
    if( test == 'gaussCItest' ) {
      indepTest<-gaussCItest
      suffStat<-list(C=cor(data),n=N,removeNAs=removeNAs)
    } else if( test == 'disCItest' ) {
      indepTest<-disCItest
      suffStat<-list(dm=data,adaptDF=FALSE)
    } else if( test == 'gaussCIcontexttest_slow' ) {
      indepTest<-gaussCIcontexttest_slow
      suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
    } else if( test == 'gaussCIcontexttest' ) {
      indepTest<-gaussCIcontexttest
      # find indices for unique joint values of context variables
      uniqueContextValues<-uniquecombs(data[,contextVars])
      regimes<-attr(uniqueContextValues,'index')
      suffStat<-list(data=data,contextVars=contextVars,uniqueContextValues=uniqueContextValues,regimes=regimes,verbose=verbose,removeNAs=removeNAs,skiptests=(mode=='jci123r'))
    } else if( test == 'gaussCIsincontest' ) {
      indepTest<-gaussCIsincontest
      suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
    } else {
      error('unknown test')
    }
  } else {
    error('unknown mode')
  }

  # run FCI
  if( mode == 'jci123' || mode == 'jci123r' || mode == 'jci1' || mode == 'jci0' ) {
    # run FCI
    if( mode == 'jci123' || mode == 'jci123r' )
      jcimode<-'123'
    else if( mode == 'jci1' )
      jcimode<-'1'
    else
      jcimode<-'0'
    pag<-fcijci(labels=colnames(data),suffStat,indepTest,alpha=alpha,verbose=verbose,doPdsep=doPDsep,contextVars=contextVars,selectionBias=FALSE,jci=jcimode)
    P<-pag@amat

    conf<-matrix(0,p,p)
    colnames(conf)<-colnames(data)
    #conf<-cpag_to_conf(P)
    
    edge<-matrix(0,p,p)
    #edge<-cpag_to_edge(P)
    if( mode == 'jci123' || mode == 'jci123r' ) {
      # direct intervention (non)-targets
      edge<-fcijci123_direct_intervention_targets(P,systemVars,contextVars)
    }
    colnames(edge)<-colnames(data)
  } else {
    pag<-fcijci(labels=colnames(data[,systemVars]),suffStat,indepTest,alpha=alpha,verbose=verbose,doPdsep=doPDsep,selectionBias=FALSE)
    # extend PAG with context variables
    P<-pag@amat
    P<-matrix(0,p,p)
    colnames(P)<-colnames(data)
    P[systemVars,systemVars]<-pag@amat
#    P[contextVars,contextVars]<-2

    # extend conf with context variables
    conf<-matrix(0,p,p)
    colnames(conf)<-colnames(data)
    #conf[systemVars,systemVars]<-cpag_to_conf(pag@amat)

    # extend edge with context variables
    edge<-matrix(0,p,p)
    colnames(edge)<-colnames(data)
    #edge[systemVars,systemVars]<-cpag_to_edge(pag@amat)
  }

  ## P: PAG
  ## P[i,j] = 0 iff no edge btw i,j
  ## P[i,j] = 1 iff i *-o j
  ## P[i,j] = 2 iff i *-> j
  ## P[i,j] = 3 iff i *-- j

  # convert PAG to ancestral relations
  arel<-cpag_to_mc(P)$C

  # convert PAG to representative MAG
#  mag<-pag2L(P)   # disabled for now because it is buggy
  mag<-list()
  mag$G<-matrix(0,p,p)
  mag$Ge<-matrix(0,p,p)
  mag$Gs<-matrix(0,p,p)
  mag$Gcircles<-matrix(0,p,p)
#  magAM <- pag2magAM(P, 1, max.chordal = 10, verbose = FALSE)  # disabled for now because it is buggy
#  mag <- magAM2mag(magAM)

  # return result
  result<-list(p=p,systemVars=systemVars,contextVars=contextVars,labels=colnames(data),pag=P,arel=arel,mag=mag,conf=conf,edge=edge)
}

fcijci123_direct_intervention_targets <- function(P,systemVars,contextVars){
  # read off direct intervention targets from PAG P produced by FCI-JCI123
  #
  # NO SELECTION BIAS ONLY! (so far)
  # ACYCLIC ONLY! (so far)
  #
  # P: PAG
  # P[i,j] = 0 iff no edge btw i,j
  # P[i,j] = 1 iff i *-o j
  # P[i,j] = 2 iff i *-> j
  # P[i,j] = 3 iff i *-- j
  
  p<-dim(P)[1]
  edge<-matrix(0,p,p)
  for( i in contextVars ) {
    for( j in systemVars ) { 
      if( P[i,j] == 0 ) {
        edge[i,j] = -1
      } else if( P[i,j] == 2 && P[j,i] == 3 ) { # i --> j
        edge[i,j] = 1
        for( k in systemVars ) {
          if( k != j ) {
            if( (P[i,k] == 2 && P[k,i] == 1) || (P[i,k] == 2 && P[k,i] == 3) || (P[i,k] == 1 && P[k,i] == 1) ) { # i o-> k, i --> k, i o-o k
              if( (P[k,j] == 2 && P[j,k] == 1) || (P[k,j] == 2 && P[j,k] == 3) || (P[k,j] == 1 && P[j,k] == 1) ) { # k o-> j, k --> j, k o-o j
                # copy P but modify k *-* j into k --> j
                Pmod<-P
                Pmod[k,j]<-2
                Pmod[j,k]<-3
                if( !visibleEdge(Pmod,k,j) ) {
                  edge[i,j] = 0
                  break
                }
              }
            }
          }
        }
      }
    }
  }

  return(edge)
}

magAM2mag<-function(magAM) {
    #magAM is a MAG (e.g., such as output pag2magAM), represented by its adjacency matrix.
    #assumes no selection bias!
    #The edge marks are encoded by numbers: 0 = no edge, 1 = circle, 2 = arrowhead, 3 = tail. If
    #amat[i,j] = 1 and amat[j,i] = 2, this represents the edge i <-o j.
  ## Notation:
  ## ----------------------------------------------------------------------
  ## amat[i,j] = 
  ##   0: no edge
  ##   1: *---o (circle)
  ##   2: *---> (arrowhead)
  ##   3: *---- (tail)
  ## (all at the j node)
   #
   # The output mag is a model as used in the HEJ2014 code
   # it has fields G, Ge, Gs, Gcircles
   #   G are the directed edges
   #   Ge the bidirected edges
   #   Gs undirected edges
   #   Gcircles are unused
   # e.g. G[to,from]=1 iff from --> to (i.e. the opposite convention as the adjacency matrix in pcalg)
    G<-magAM
    N<-dim(G)[1]
    mag<-list(G=array(0,c(N,N)),
              Ge=array(0,c(N,N)),
              Gs=array(0,c(N,N)),
              Gcircles=array(0,c(N,N)))
    #now must turn G into G, Ge, Gs, Gcircles matrix
    #note that the rows and cols in the two notations are switched
    
    while (any(G!=0) ) {
      vars<-which(G!=0,arr.ind=TRUE)[1,]
      from<-vars[1];to<-vars[2];

      if( G[to,from] == 1 || G[from,to] == 1 ) {
        stop('magAM must not have circles in it')
      } else if( G[to,from] == 3 && G[from,to] == 3 ) {#type == '---'
        stop('undirected edge in magAM')
      } else if( G[to,from] == 2 &&  G[from,to] == 2 ) {#type == '<->'
        mag$Ge[to,from]<-mag$Ge[from,to]<-1
      } else if( G[to,from] == 3 && G[from,to] == 2 ) {#type == '-->'
        mag$G[to,from]<-1
      } else if( G[to,from] == 2 && G[from,to] == 3 ) {#type == '<--'
        mag$G[from,to]<-1
      } else {
        stop('unknown edge type in magAM')
      }
      
      G[from,to]<-G[to,from]<-0
    }
  return(mag)
}
