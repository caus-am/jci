# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

icp_wrapper <- function(data,systemVars,contextVars,alpha=1e-2,verbose=0,subsamplefrac=0.0,datamode='merge') {
  # data:          Nxp matrix
  # systemVars:    indices of sysem variables (in 1..p)
  # contextVars:   indices of context variables (in 1..p)
  # alpha:         p-value threshold for independence test (default: 1e-2)
  # verbose:       verbosity level (default: 0)
  # subsamplefrac: fraction of subsamples used for subsampling; if 0, don't subsample and use all data
  # datamode:      one of {'multiple','merge'}, default='merge'

  suppressMessages(library(InvariantCausalPrediction))

  # prepare data
  X<-prepare_jci(data=data,systemVars=systemVars,contextVars=contextVars,subsamplefrac=subsamplefrac,mode=datamode,obsContext=c())
  data<-X$data
  systemVars<-X$systemVars
  contextVars<-X$contextVars
  N<-X$N
  p<-X$p
  removeNAs<-X$removeNAs

  pSys <- length(systemVars)
  pCon <- length(contextVars)
  arel <- matrix(0,p,p)
  edge <- matrix(0,p,p)
  conf <- matrix(0,p,p)
  colnames(arel) <- colnames(data)
  colnames(edge) <- colnames(data)
  colnames(conf) <- colnames(data)

  if( pCon > 0 ) {
    # run ICP
    if( datamode == 'merge' ) {
      for( target in 1:pSys ) {
        if( verbose >= 2 )
          cat('Target variable: ',target,'\n')
        result <- ICP(as.matrix(data[,setdiff(systemVars,target)]),data[,target],ExpInd=data[,contextVars],alpha=alpha,showCompletion=FALSE,showAcceptedSets=verbose>=2,stopIfEmpty=TRUE)
        if( verbose >= 2 )
          cat('\n')
        # look at results only if not whole model has been rejected
        if( !result$modelReject ) {
          # check if significant results
          if( any(result$pvalues <= alpha) ) {
            causes<-which(result$pvalues<=alpha)
            idx<-setdiff(1:pSys,target)
            arel[idx[causes], target] <- -log(result$pvalues[causes])
          }
        } else {
          if( verbose >= 2 ) {
            cat('Model rejected\n')
            cat(result$pvalues,'\n')
          }
        }
      }
    } else if( datamode == 'multiple' ) {
      for( target in 1:pSys ) {
        for( cV in 1:pCon ) {
          if( verbose >= 2 )
            cat('Target variable: ',target,', context variable: ',contextVars[cV],'\n')
          result <- ICP(as.matrix(data[,setdiff(systemVars,target)]),data[,target],ExpInd=data[,contextVars[cV]],alpha=alpha,showCompletion=FALSE,showAcceptedSets=verbose>=2,stopIfEmpty=TRUE)
          if( verbose >= 2 )
            cat('\n')
          # look at results only if not whole model has been rejected
          if( !result$modelReject ) {
            # check if significant results
            if( any(result$pvalues <= alpha) ) {
              causes<-which(result$pvalues<=alpha)
              idx<-setdiff(1:pSys,target)
              #arel[causes, target] <- -log(result$pvalues[causes])
              for( c in causes ) {
                arel[idx[c],target] <- max(arel[idx[c],target],-log(result$pvalues[c]))
              }
            }
          } else {
            if( verbose >= 2 ) {
              cat('Model rejected\n')
              cat(result$pvalues,'\n')
            }
          }
        }
      }
    }
  }

  # return result
  # ICP returns ancestral relations in the causally insufficient case, and direct relations in the causally sufficient case.
  result<-list(p=p,systemVars=systemVars,contextVars=contextVars,labels=colnames(data),arel=arel,edge=edge,conf=conf)
}
