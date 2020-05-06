# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

asd_wrapper <- function(indepfile,data,systemVars,contextVars,alpha=1e-2,verbose=0,subsamplefrac=0.0,mode='obs',test='gaussCIcontexttest',obsContext=matrix(0,1,length(contextVars)),weightmul=1000,extrafiles,Itargets=c()) {
  # indepfile:     name of ASP file with independence test results
  # data:          Nxp matrix
  # systemVars:    indices of sysem variables (in 1..p)
  # contextVars:   indices of context variables (in 1..p)
  # alpha:         p-value threshold for independence test (default: 1e-2)
  # verbose:       verbosity level (default: 0)
  # subsamplefrac: fraction of subsamples used for subsampling; if 0, don't subsample and use all data
  # mode:          'obs' / 'pooled' / 'meta' / 'jci13' / 'jci1' / 'jci1nt' / 'jci12' / 'jci12nt' / 'jci123kt' / 'pikt' / 'jci123-sc' / 'jci1-sc' / 'jci0' / 'jci0nt' / 'jci123io'
  # test:          'gaussCItest' / 'gaussCIcontexttest_slow' / 'gaussCIcontexttest' / 'gaussCIsincontest' (only relevant for modes 'jci123','jci13','jci1','jci1nt','jci123kt','jci12','jci12nt','jci0','jci0nt','jci123io')
  #                   gaussCItest: just partial correlations
  #                   gaussCIcontexttest_slow: tests within each context value separately and combines with Fisher's method
  #                   gaussCIcontexttest: tests within each context value separately and combines with Fisher's method, faster implementation
  # obsContext:    vector of values of context variables that defines the observational context
  # weightmul:     weight multiplier for ASP (default: 1000)
  # extrafiles:    additional ASP files
  # Itargets:      intervention targets of perfect interventions for each regime (if asd-pikt needs to be run)

  suppressMessages(library(pcalg))
  suppressMessages(library(mgcv))

  # prepare data
  if( mode == 'obs' )
    X<-prepare_jci(data,systemVars,contextVars,subsamplefrac,'obsonly',obsContext)
  else if( mode == 'jci123-sc' || mode == 'jci1-sc' )
    X<-prepare_jci(data,systemVars,contextVars,subsamplefrac,'merge',obsContext=c())
  else
    X<-prepare_jci(data,systemVars,contextVars,subsamplefrac,'multiple',obsContext=c())
  data<-X$data
  systemVars<-X$systemVars
  contextVars<-X$contextVars
  N<-X$N
  p<-X$p
  removeNAs<-X$removeNAs
  labels<-colnames(data)

  # setup independence test
  if( mode == 'obs' || mode == 'pooled' || length(contextVars) == 0 ) {
    if( removeNAs )
      stop('removeNAs not implemented yet for this mode')
    indepTest<-gaussCItest
    if( mode == 'pikt' ) {
      Cs<-vector(mode='list',length=1)
      ns<-vector(mode='list',length=1)
      Cs[[1]]<-cor(data[,systemVars])
      ns[[1]]<-N
      suffStat<-list(Cs=Cs,ns=ns)
    } else
      suffStat<-list(C=cor(data[,systemVars]),n=N)
  } else if( mode == 'meta' || mode == 'pikt' ) {
    if( removeNAs )
      stop('removeNAs not implemented yet for this mode')
    # find indices for unique joint values of context variables
    uniqueContextValues<-uniquecombs(data[,contextVars])
    regimes<-attr(uniqueContextValues,'index')
    nRegimes<-max(regimes)

    if( mode == 'meta' )
      indepTest<-gaussCIFishertest # use Fisher's method for combining p-values tested in different regimes
    else if( mode == 'pikt' )
      indepTest<-gaussCItest
    Cs<-vector(mode='list',length=nRegimes)
    ns<-vector(mode='integer',length=nRegimes)
    for( R in 1:nRegimes ) {
      subset<-which(regimes==R)
      Cs[[R]]<-cor(data[subset,systemVars])
      ns[[R]]<-length(subset)
    }
    suffStat<-list(Cs=Cs,ns=ns,verbose=verbose)
  } else if( mode == 'jci123-sc' || mode == 'jci1-sc' ) {
    indepTest<-gaussCIsincontest
    suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
  } else if( mode == 'jci123' || mode == 'jci13' || mode == 'jci1' || mode == 'jci1nt' || mode == 'jci123kt' || mode == 'jci12' || mode == 'jci12nt' || mode == 'jci0' || mode == 'jci0nt' || mode == 'jci123io' ) {
    if( test == 'gaussCItest' ) {
      indepTest<-gaussCItest
      suffStat<-list(C=cor(data),n=N,removeNAs=removeNAs)
    } else if( test == 'gaussCIcontexttest_slow' ) {
      indepTest<-gaussCIcontexttest_slow
      suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
    } else if( test == 'gaussCIsincontest' ) {
      indepTest<-gaussCIsincontest
      suffStat<-list(data=data,contextVars=contextVars,verbose=verbose,removeNAs=removeNAs)
    } else if( test == 'gaussCIcontexttest' ) {
      indepTest<-gaussCIcontexttest
      # find indices for unique joint values of context variables
      uniqueContextValues<-uniquecombs(data[,contextVars])
      regimes<-attr(uniqueContextValues,'index')
      suffStat<-list(data=data,contextVars=contextVars,uniqueContextValues=uniqueContextValues,regimes=regimes,verbose=verbose,removeNAs=removeNAs)
    } else {
      stop('unknown test')
    }
  } else {
    stop('unknown mode')
  }

  # run ASD
  if( mode == 'jci123' || mode == 'jci13' || mode == 'jci1' || mode == 'jci1nt' || mode == 'jci123kt' || mode == 'jci12' || mode == 'jci12nt' || mode == 'jci123-sc' || mode == 'jci1-sc' || mode == 'jci0' || mode == 'jci0nt' || mode == 'jci123io' ) {
    # JCI modes
    asd_indeptests(indepfile,data,contextVars,alpha=alpha,verbose=verbose,indepTest,suffStat,weightmul,alltests=(mode=='jci1') || (mode=='jci12') || (mode=='jci0'),iotests=(mode == 'jci123io'))

    if( mode == 'jci123' || mode == 'jci123kt' || mode == 'jci123-sc' || mode == 'jci123io' ) 
      extrafiles<-c(extrafiles,file.path(asp_dir,'jci123.pl'))
    else if( mode == 'jci13' )
      extrafiles<-c(extrafiles,file.path(asp_dir,'jci1.pl'),file.path(asp_dir,'jci3.pl'))
    else if( mode == 'jci1' || mode == 'jci1nt' || mode == 'jci1-sc' )
      extrafiles<-c(extrafiles,file.path(asp_dir,'jci1.pl'))
    else if( mode == 'jci12' || mode == 'jci12nt' )
      extrafiles<-c(extrafiles,file.path(asp_dir,'jci1.pl'),file.path(asp_dir,'jci2.pl'))
    #else if( mode == 'jci0' || mode == 'jci0nt' )
      #nothing to add

    # query ancestors, edges and confounders
    result<-list(p=p,systemVars=systemVars,contextVars=contextVars,labels=labels)
    result$arel<-query_clingo_arel(p,fromVars=1:p,toVars=1:p,basefilename=indepfile,extrafiles=extrafiles,labels=labels)
    result$edge<-query_clingo_edge(p,fromVars=1:p,toVars=1:p,basefilename=indepfile,extrafiles=extrafiles,labels=labels)
    result$conf<-query_clingo_conf(p,fromVars=1:p,toVars=1:p,basefilename=indepfile,extrafiles=extrafiles,labels=labels)
  } else { # mode: 'obs' / 'pooled' / 'meta' / 'pikt'
    # nonJCI modes
    if( mode == 'pikt' ) {
      warning('hacky implementation!')
      asd_pikt_indeptests(indepfile,data[,systemVars],contextVars=c(),alpha=alpha,verbose=verbose,indepTest,suffStat,weightmul,Itargets)
    } else {
      asd_indeptests(indepfile,data[,systemVars],contextVars=c(),alpha=alpha,verbose=verbose,indepTest,suffStat,weightmul,alltests=FALSE,iotests=FALSE)
    }

    # query ancestors, edges and confounders
    pSys<-length(systemVars)
    result<-list(p=p,systemVars=systemVars,contextVars=contextVars,labels=labels)
    result$arel<-query_clingo_arel(p,fromVars=1:pSys,toVars=1:pSys,basefilename=indepfile,extrafiles=extrafiles,labels=labels)
    result$edge<-query_clingo_edge(p,fromVars=1:pSys,toVars=1:pSys,basefilename=indepfile,extrafiles=extrafiles,labels=labels)
    result$conf<-query_clingo_conf(p,fromVars=1:pSys,toVars=1:pSys,basefilename=indepfile,extrafiles=extrafiles,labels=labels)
  }

  # return result
  return(result)
}
