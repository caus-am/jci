# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

simul_linear_stochastic <- function(basefilename,pSysObs,pContext,eps,eta,N,acyclic,surgical,seed=0) {
# simulate from linear-Gaussian SCM with diagonal design, stochastic interventions
#
# order variables as: observed system variables, context variables, latent confounding system variables

suppressMessages(library(MASS))
suppressMessages(library(rjson))
suppressMessages(library(expm))
source('../../init.R',chdir=TRUE)

stopifnot( pSysObs > 0 )

# set random number seed
set.seed(seed)

# choose bidirected edges randomly with density eta
U <- matrix(as.numeric(runif(pSysObs*pSysObs)<eta),pSysObs,pSysObs)
U <- U - lower.tri(U,diag=TRUE) * U
pSysConf <- length(which(U!=0))

# complete adjacency matrix
p <- pSysObs + pContext + pSysConf
A <- matrix(0,p,p)

# choose causal relations between observed system variables randomly with density eps
repeat {
  A[1:pSysObs,1:pSysObs] <- matrix(as.numeric(runif(pSysObs*pSysObs)<eps),pSysObs,pSysObs)
  # remove (self)loops
  if( acyclic ) {
    A <- A - lower.tri(A,diag=TRUE) * A
    Aacyclic <- TRUE
    break
  } else {
    A <- A - diag(diag(A))
    # test for acyclicity
    Aacyclic <- sum((diag(expm(A[1:pSysObs,1:pSysObs]))-rep(1,pSysObs))^2) < 1e-10
    if( !Aacyclic )
      break
  }
}

# use one confounder for each bidirected edge
if( pSysConf > 0 ) {
  Utargets <- which(U!=0,arr.ind=T)
  for( i in 1:pSysConf ) {
    A[pSysObs + pContext + i,Utargets[i,]] <- 1
  }
}

# choose one intervention target for each context variable
# if pContext <= pSysObs, ensure that no target is shared by multiple context variables
targets <- sample(pSysObs,size=pContext,replace=(pContext > pSysObs))
if( pContext > 0 )
  for( i in 1:pContext )
    A[pSysObs + i,targets[i]] <- 1

# choose weights uniformly in [-1.5,-0.5] \cup [0.5,1.5]
# weights <- matrix(rnorm(p*p),p,p)
weights <- sign(matrix(rnorm(p*p),p,p)) * (matrix(runif(p*p),p,p) + 0.5)
B <- A * weights

# rescale weight matrix so that every system variable has more or less similar scale
# we rescale such that if all inputs would be i.i.d. standard-normal, the output also has variance 1
# this is an approximation because it ignores covariances between inputs
scales<-matrix(0,p,1)
if( pSysConf == 0 )
  sysVars<-c(1:pSysObs)
else
  sysVars<-c(1:pSysObs,(pSysObs+pContext+1):p)
for( i in 1:pSysObs ) {
  scales[i] <- (sum(B[sysVars,i]^2) + 1)
  B[,i] <- B[,i] / scales[i]^0.5
}

# write weight matrix
write.csv(B,file=paste(basefilename,'-B.csv',sep=''),row.names=FALSE)

# save variable types
systemVars <- 1:pSysObs
if( pContext > 0 )
  contextVars <- (pSysObs+1):(pSysObs+pContext)
else
  contextVars <- c()

# make list of regimes ('diagonal' design)
nRegimes <- 1 + pContext   # also an observational one
# experimental design
expdesign<-matrix(0,nRegimes,pContext)
if( nRegimes >= 2 ) 
  for( R in 2:nRegimes )
    expdesign[R,R-1] <- 1

# simulate data for each regime
Ns <- vector("list",nRegimes)
Xs <- vector("list",nRegimes)
Es <- vector("list",nRegimes)
Bs <- vector("list",nRegimes)
Itargets <- vector("list",nRegimes)
for( R in 1:nRegimes ) {
  # number of data points to simulate
  Ns[[R]] <- N

  # noise variables
#  Es[[R]] <- matrix(runif(N*p),N,p) - 0.5
#  Es[[R]] <- matrix(rlaplace(N*p),N,p)
  Es[[R]] <- matrix(rnorm(Ns[[R]]*p),Ns[[R]],p)
  # override context noises with exogenous choice for this regime
  if( pContext > 0 )
    Es[[R]][,contextVars] <- matrix(expdesign[R,],Ns[[R]],length(contextVars),byrow=T)

  # find union of intervened variables
  Iset <- c()
  if( pContext > 0 )
    for( i in 1:pContext )
      if( expdesign[R,i] )
        Iset <- union(Iset, which(A[pSysObs + i,] != 0))
  # save intervention targets as integer in binary encoding
  Itargets[[R]] <- 0
  for( b in 1:p )
    Itargets[[R]] <- (b %in% Iset)*2^(b-1) + Itargets[[R]]

  # choose values for intervened variables
#  Ivals <- matrix(runif(N*length(Iset)),N,length(Iset)) - 0.5
#  Ivals <- matrix(rlaplace(N*length(Iset)),N,length(Iset))
  Ivals <- matrix(rnorm(Ns[[R]]*length(Iset)),Ns[[R]],length(Iset)) + 1
  Es[[R]][,Iset] <- Ivals

  # copy B and cut incoming edges for intervened variables
  Bs[[R]] <- B
  if( surgical )
    Bs[[R]][,Iset] <- 0

  # solve for X
  Xs[[R]] <- t(solve(diag(p)-t(Bs[[R]]), t(Es[[R]])))
}

# concatenate everything
X <- matrix(0,0,p)
for( R in 1:nRegimes )
  X <- rbind(X,Xs[[R]])

# write results to file
write.csv(X,file=paste(basefilename,'-data.csv',sep=''),row.names=FALSE)

# save augmented causal graph
write.csv(A,file=paste(basefilename,'-graphaug.csv',sep=''),row.names=FALSE)
dg2graphviz(filename=paste(basefilename,'-graphaug.dot',sep=''),p=p,contextVars=contextVars,labels=1:p,dg=A)

# construct marginal causal graph on observed variables in HEJ2014 convention
L <- list()
L$G <- t(A[1:(pSysObs+pContext),1:(pSysObs+pContext)])
L$Ge <- matrix(0,pSysObs+pContext,pSysObs+pContext)
if( pSysConf > 0 )
  for( i in 1:pSysConf )
    for( j1 in 1:pSysObs )
      if( A[pSysObs+pContext+i,j1] )
        if( j1 < pSysObs ) {
          for( j2 in (j1+1):pSysObs )
            if( A[pSysObs+pContext+i,j2] ) {
              L$Ge[j1,j2] <- 1
              L$Ge[j2,j1] <- 1
            }
        }
L$Gs <- matrix(0,pSysObs+pContext,pSysObs+pContext)
L$Gcircles <- matrix(0,pSysObs+pContext,pSysObs+pContext)
L2graphviz(filename=paste(basefilename,'-graph.dot',sep=''),p=pSysObs+pContext,contextVars=contextVars,labels=1:(pSysObs+pContext),L=L)

# save edges
edge<-A[1:(pSysObs+pContext),1:(pSysObs+pContext)]
write.csv(edge,file=paste(basefilename,'-edge.csv',sep=''),row.names=FALSE)

# save confounders
write.csv(L$Ge,file=paste(basefilename,'-conf.csv',sep=''),row.names=FALSE)

# save ancestral relations
anc<-sign(expm(A[1:(pSysObs+pContext),1:(pSysObs+pContext)]))
anc<-anc-diag(diag(anc))
write.csv(anc,file=paste(basefilename,'-arel.csv',sep=''),row.names=FALSE)

# write metadata
metadata<-list()
metadata$SystemVars<-as.list(systemVars)
metadata$ContextVars<-as.list(contextVars)
metadata$basefilename<-basefilename
metadata$pSysObs<-pSysObs
metadata$pSysConf<-pSysConf
metadata$pContext<-pContext
metadata$eps<-eps
metadata$eta<-eta
metadata$N<-N*nRegimes
metadata$acyclic<-acyclic
metadata$sufficient<-0
metadata$surgical<-surgical
metadata$seed<-seed
metadata$Itargets<-Itargets
metadata$targets<-as.list(2^(targets-1))

cat(file=paste(basefilename,'.json',sep=''),toJSON(metadata))
}
