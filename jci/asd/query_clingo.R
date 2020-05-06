# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

query_clingo_arel <- function(p,fromVars,toVars,basefilename,extrafiles,clingo_cmd='clingo',clingo_pars='--quiet=2,1 -W no-atom-undefined',verbose=0,labels=c()) {

  arel<-matrix(0,p,p)
  stopifnot(min(fromVars) >= 1 && max(fromVars) <= p)
  stopifnot(min(toVars) >= 1 && max(toVars) <= p)
  for( i in fromVars )
    for( j in toVars ) 
      if( i != j ) {
        val0 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=paste('ancestor(',i-1,',',j-1,',0).',sep=''))
        val1 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=paste(':-ancestor(',i-1,',',j-1,',0).',sep=''))
        arel[i,j] <- val1 - val0
      }

  if( !is.null(labels) )
    colnames(arel)=labels
  return(arel)
}

query_clingo_edge <- function(p,fromVars,toVars,basefilename,extrafiles,clingo_cmd='clingo',clingo_pars='--quiet=2,1 -W no-atom-undefined',verbose=0,labels=c()) {

  edge<-matrix(0,p,p)
  stopifnot(min(fromVars) >= 1 && max(fromVars) <= p)
  stopifnot(min(toVars) >= 1 && max(toVars) <= p)
  for( i in fromVars )
    for( j in toVars ) 
      if( i != j ) {
        val0 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=paste('edge(',i-1,',',j-1,').',sep=''))
        val1 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=paste(':-edge(',i-1,',',j-1,').',sep=''))
        edge[i,j] <- val1 - val0
      }

  if( !is.null(labels) )
    colnames(edge)=labels
  return(edge)
}

query_clingo_conf <- function(p,fromVars,toVars,basefilename,extrafiles,clingo_cmd='clingo',clingo_pars='--quiet=2,1 -W no-atom-undefined',verbose=0,labels=c()) {

  conf<-matrix(0,p,p)
  stopifnot(min(fromVars) >= 1 && max(fromVars) <= p)
  stopifnot(min(toVars) >= 1 && max(toVars) <= p)
  for( i in fromVars )
    for( j in toVars ) 
      if( i < j ) {
        val0 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=paste('conf(',i-1,',',j-1,').',sep=''))
        val1 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=paste(':-conf(',i-1,',',j-1,').',sep=''))
        conf[i,j] <- val1 - val0
        conf[j,i] <- conf[i,j]
      }

  if( !is.null(labels) )
    colnames(conf)=labels
  return(conf)
}

query_clingo_indep <- function(p,fromVars,toVars,basefilename,extrafiles,clingo_cmd='clingo',clingo_pars='--quiet=2,1 -W no-atom-undefined',verbose=0,labels=c()) {
#  conf<-matrix(0,p,p)
  stopifnot(min(fromVars) >= 1 && max(fromVars) <= p)
  stopifnot(min(toVars) >= 1 && max(toVars) <= p)
  result=''
  for( i in fromVars )
    for( j in toVars ) 
      if( i < j ) {
        if( i != 1 || j != 2 ) {
          for( cset in 0:(2^p-1) ) {
            if( !bitwAnd(cset,2^(i-1)) && !bitwAnd(cset,2^(j-1)) ) {
              kset<-(2^p-1) - 2^(i-1) - 2^(j-1) - cset
              strij<-paste(i-1,',',j-1,',',cset,',0,',kset,sep='')
              strji<-paste(j-1,',',i-1,',',cset,',0,',kset,sep='')
              # add dep statement
              query0 <- paste(':- not th(',strij,'), not th(',strji,'), not hh(',strij,'), not tt(',strij,').',sep='')
              val0 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=query0)
              # add indep statement
              query1 <- paste(':- th(',strij,'). :- th(',strji,'). :- hh(',strij,'). :- tt(',strij,').',sep='')
              val1 <- run_clingo(basefilename=basefilename,extrafiles=extrafiles,clingo_cmd=clingo_cmd,clingo_pars=clingo_pars,verbose=verbose,query=query1)
              if( verbose == 2 ) 
                cat(i-1,' indep ',j-1,' given ',cset,'? ',val1 - val0,'\n',sep='')
              if( val0 == Inf ) {
                if( val1 == Inf ) {
                  result=paste(result,'% CONTRADICTION: dep(',strij,',Inf) AND indep(',strij,',Inf).\n',sep='')
                } else {
                  #result=paste(result,'% indep(',strij,',Inf).\n:- th(',strij,'). :- th(',strji,'). :- hh(',strij,'). :- tt(',strij,').\n',sep='')
                  result=paste(result,'indep(',strij,',1).\n',sep='')
                }
              } else {
                if( val1 == Inf ) {
#                  result=paste(result,'% dep(',strij,',Inf).\n:- not th(',strij,'), not th(',strji,'), not hh(',strij,'), not tt(',strij,').\n',sep='')
                  result=paste(result,'dep(',strij,',1000).\n',sep='')
                } else {
                  if( val0 > val1 ) {
                    result=paste(result,'indep(',strij,',',val0 - val1,').\n',sep='')
                  } else if( val0 < val1 ) {
                    result=paste(result,'dep(',strij,',',val1 - val0,').\n',sep='')
                  } else {
                    result=paste(result,'%??dep(',strij,').\n',sep='')
                  }
                }
              }
      #        conf[i,j] <- val1 - val0
      #        conf[j,i] <- conf[i,j]
            }
          }
        }
      }

#  if( !is.null(labels) )
#    colnames(conf)=labels
#  return(conf)
  return(result)
}
