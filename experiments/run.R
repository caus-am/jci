# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# Use this script to run a method on simulated data

suppressMessages(library(rjson))
source('../../init.R',chdir=TRUE)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

basefilename<-args[1]
# read metadata
metadatafile=paste(basefilename,'.json',sep='')
metadata<-fromJSON(file=metadatafile)
#pSysObs<-metadata$pSysObs
#stopifnot( pSysObs > 0 )
#pContext<-metadata$pContext
systemVars<-metadata$SystemVars
contextVars<-metadata$ContextVars

pSysObs<-length(systemVars)
stopifnot( pSysObs > 0 )
pContext<-length(contextVars)
p <- pSysObs + pContext
obsContext <- matrix(0,1,pContext)

#alg<-'fci'
alg<-args[2]
#mode<-'obs'
mode<-args[3]
miniter<-0
if( !is.na(args[4]) ) {
  miniter<-as.numeric(args[4])
}
maxiter<-0
if( !is.na(args[5]) ) {
  maxiter<-as.numeric(args[5])
}
alpha<-1e-2
if( !is.na(args[6]) )
  alpha<-as.numeric(args[6])
jcifci_test<-'gaussCIcontexttest'
if( !is.na(args[7]) )
  jcifci_test<-args[7]
doPDsep<-TRUE
if( !is.na(args[8]) )
  doPDsep<-(as.numeric(args[8]) != 0)
verbose<-0
if( !is.na(args[9]) )
  verbose<-(as.integer(args[9]))

cat('run.R: running with arguments',args,'\n')

# read and preprocess data
data<-read.csv(file=paste(basefilename,'-data.csv',sep=''),header=TRUE,sep=",")

for( iter in miniter:maxiter ) {
  cat('iter: ',iter,'\n')

set.seed(iter)
if( iter != 0 ) {
  subsamplefrac<-0.5
} else {
  subsamplefrac<-0.0
}

# start measuring time
start_time<-proc.time()[3]

# run causal-discovery
if( alg=='asd' ) { # run ASD
  if( (mode == 'obs' && pSysObs <= 6) || 
      (mode == 'pooled' && pSysObs <= 6) || 
      (mode == 'meta' && pSysObs <= 6) || 
      (mode == 'pikt' && pSysObs <= 6) || 
      (mode == 'jci123' && p <= 8) || 
      (mode == 'jci123kt' && p <= 8) ||
      (mode == 'jci13' && p <= 6) ||
      (mode == 'jci1' && p <= 6) ||
      (mode == 'jci1nt' && p <= 6) ||
      (mode == 'jci12' && p <= 6) ||
      (mode == 'jci12nt' && p <= 6) ||
      (mode == 'jci123-sc' && p <= 7) ||
      (mode == 'jci1-sc' && p <= 7) ||
      (mode == 'jci0' && p <= 6) ||
      (mode == 'jci0nt' && p <= 6) ||
      (mode == 'jci123io' && p <= 8)
     ) {
    if( iter != 0 )
      outfile<-paste(basefilename,'-',alg,'-',mode,'-',iter,sep='')
    else
      outfile<-paste(basefilename,'-',alg,'-',mode,sep='')
    indepfile<-paste(outfile,'.indep',sep='')

    if( mode=='pikt' ) # reason about perfect interventions
      extrafiles=c(file.path(asp_dir,'partial_comp_tree.pl'))
    else # can forget about perfect interventions
      extrafiles=c(file.path(asp_dir,'obs_comp_tree.pl'))
    if( metadata$acyclic ) { # use acyclic d-separation encoding
      if( is.null(metadata$sufficient) || !metadata$sufficient )
        extrafiles=c(extrafiles,file.path(asp_dir,'asd_acyclic.pl'))
      else
        extrafiles=c(extrafiles,file.path(asp_dir,'asd_acyclic_sufficient.pl'))
    } else { # cyclic, so run with sigma-separation (note: the model as a whole is not linear)
      extrafiles=c(extrafiles,file.path(asp_dir,'asd_sigma_cyclic.pl'))
    }

    if( mode=='jci123kt' ) {
      targetfile<-paste(outfile,'.targets',sep='')
      f=file(targetfile,'w')
      if( length(contextVars) > 0 )
        for( i in 1:length(contextVars) ) 
          for( j in 1:length(systemVars) ) {
            targets <- which(intToBits(metadata$targets[i])!=0)
            if( j %in% targets ) {
              cat(file=f,'edge(',contextVars[i]-1,',',systemVars[j]-1,').\n')
            } else {
              cat(file=f,':-edge(',contextVars[i]-1,',',systemVars[j]-1,').\n')
            }
          }
      close(f)
      extrafiles=c(extrafiles,targetfile)
    }

    result<-asd_wrapper(indepfile,data,systemVars,contextVars,alpha,verbose=verbose,subsamplefrac,mode,test='gaussCIcontexttest',obsContext,weightmul=1000,extrafiles,Itargets=metadata$Itargets)
    write.csv(result$arel,file=paste(outfile,'-arel.csv',sep=''),row.names=FALSE)
    write.csv(result$edge,file=paste(outfile,'-edge.csv',sep=''),row.names=FALSE)
    write.csv(result$conf,file=paste(outfile,'-conf.csv',sep=''),row.names=FALSE)
  } else
    cat('Skipping asd-',mode,' because unknown mode or p=',p,'\n',sep='')
} else if( alg == 'fci' ) { # run FCI
  result<-fci_wrapper(data,systemVars,contextVars,alpha,verbose=verbose,subsamplefrac,test=jcifci_test,mode,obsContext,doPDsep)

  # save results
  #ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)
  if( iter != 0 )
    outfile<-paste(basefilename,'-',alg,'-',mode,'-',iter,sep='')
  else
    outfile<-paste(basefilename,'-',alg,'-',mode,sep='')
  save(result,file=paste(outfile,'-result.Rdata',sep=''))
  pag2graphviz(filename=paste(outfile,'-pag.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,pag=result$pag,mode=mode)
  mc2graphviz(filename=paste(outfile,'-arel.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,mc=result$arel)
  L2graphviz(filename=paste(outfile,'-mag.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,L=result$mag,mode=mode)
  write.csv(result$arel,file=paste(outfile,'-arel.csv',sep=''),row.names=FALSE)
  write.csv(result$edge,file=paste(outfile,'-edge.csv',sep=''),row.names=FALSE)
  write.csv(result$conf,file=paste(outfile,'-conf.csv',sep=''),row.names=FALSE)

  # demo: query an independence relation
  # cat('1 _||_ 11 | 2? ', !directed_reachable(1,11,c(2),c(),result$mag,verbose=0),'\n')
} else if( alg == 'icp' ) { # run ICP
  if( mode == 'mc' )
    datamode<-'multiple'
  else if( mode == 'sc' )
    datamode<-'merge'
  result<-icp_wrapper(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac=subsamplefrac,datamode=datamode)

  if( iter != 0 )
    outfile<-paste(basefilename,'-',alg,'-',mode,'-',iter,sep='')
  else
    outfile<-paste(basefilename,'-',alg,'-',mode,sep='')
  write.csv(result$arel,file=paste(outfile,'-arel.csv',sep=''),row.names=FALSE)
  save(result,file=paste(outfile,'-result.Rdata',sep=''))
  mc2graphviz(filename=paste(outfile,'-arel.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,mc=result$arel)
} else if( alg == 'lcd' ) { # run LCD
  if( mode == 'mc' ) 
    result<-lcd(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCItest',conservative=FALSE)
  else if( mode == 'mccon' ) 
    result<-lcd(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCItest',conservative=TRUE)
  else if( mode == 'mcsct' )
    result<-lcd(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCIsincontest',conservative=FALSE)
  else if( mode == 'mcsctcon' )
    result<-lcd(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCIsincontest',conservative=TRUE)
  else if( mode == 'sc' )
    result<-lcd(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='merge',test='gaussCIsincontest',conservative=FALSE)
  else if( mode == 'sccon' )
    result<-lcd(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='merge',test='gaussCIsincontest',conservative=TRUE)
  else
    stop('unknown lcd mode')

  if( iter != 0 )
    outfile<-paste(basefilename,'-',alg,'-',mode,'-',iter,sep='')
  else
    outfile<-paste(basefilename,'-',alg,'-',mode,sep='')
  save(result,file=paste(outfile,'-result.Rdata',sep=''))
  write.csv(result$arel,file=paste(outfile,'-arel.csv',sep=''),row.names=FALSE)
  mc2graphviz(filename=paste(outfile,'-arel.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,mc=result$arel > 0)
  write.csv(result$conf,file=paste(outfile,'-conf.csv',sep=''),row.names=FALSE)
  write.csv(result$edge,file=paste(outfile,'-edge.csv',sep=''),row.names=FALSE)
} else if( alg == 'cif' ) { # run CIF
  if( mode == 'mc' ) 
    result<-cif(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCItest',conservative=FALSE,patterns=c(0,1))
  else if( mode == 'mccon' ) 
    result<-cif(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCItest',conservative=TRUE,patterns=c(0,1))
  else if( mode == 'mcsct' ) 
    result<-cif(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCIsincontest',conservative=FALSE,patterns=c(0,1))
  else if( mode == 'mcsctcon' ) 
    result<-cif(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='multiple',test='gaussCIsincontest',conservative=TRUE,patterns=c(0,1))
  else if( mode == 'sc' )
    result<-cif(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='merge',test='gaussCIsincontest',conservative=FALSE,patterns=c(0,1))
  else if( mode == 'sccon' )
    result<-cif(data,systemVars,contextVars,alpha=alpha,verbose=verbose,subsamplefrac,datamode='merge',test='gaussCIsincontest',conservative=TRUE,patterns=c(0,1))
  else
    stop('unknown cif mode')

  if( iter != 0 )
    outfile<-paste(basefilename,'-',alg,'-',mode,'-',iter,sep='')
  else
    outfile<-paste(basefilename,'-',alg,'-',mode,sep='')
  save(result,file=paste(outfile,'-result.Rdata',sep=''))
  write.csv(result$arel,file=paste(outfile,'-arel.csv',sep=''),row.names=FALSE)
  mc2graphviz(filename=paste(outfile,'-arel.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,mc=result$arel)
  write.csv(result$conf,file=paste(outfile,'-conf.csv',sep=''),row.names=FALSE)
} else if( alg == 'fisher' ) { # run Fisher's method
  result<-fisher(data,systemVars,contextVars,alpha,verbose=verbose,subsamplefrac)

  if( iter != 0 )
    outfile<-paste(basefilename,'-',alg,'-',iter,sep='')
  else
    outfile<-paste(basefilename,'-',alg,sep='')
  save(result,file=paste(outfile,'-result.Rdata',sep=''))
  write.csv(result$arel,file=paste(outfile,'-arel.csv',sep=''),row.names=FALSE)
  mc2graphviz(filename=paste(outfile,'-arel.dot',sep=''),p=result$p,contextVars=result$contextVars,labels=result$labels,mc=result$arel > 0)
} else
  stop('Unknown algorithm')

# stop measuring time
stop_time<-proc.time()[3]

if( exists('outfile') ) {
  # write time to file
  cat(file=paste(outfile,'.runtime',sep=''),stop_time-start_time,'\n')
}

}
