# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# Use this script to gather the results of the bootstrap procedure

analyze <- function(basefilename,algmodes,type) {
  suppressMessages(library(rjson))
  suppressMessages(library(tools))
  source('../../init.R',chdir=TRUE)

  stopifnot(type == 'arel' || type == 'edge' || type == 'conf' || type == 'runtime' || type == 'cleanup')

  metadatafile=paste(basefilename,'.json',sep='')
  metadata<-fromJSON(file=metadatafile)
  systemVars<-metadata$SystemVars
  stopifnot( length(systemVars) > 0 )
  contextVars<-metadata$ContextVars
  obsContext <- metadata$obsContext

  pSysObs<-length(systemVars)
  pContext<-length(contextVars)
  p <- pSysObs + pContext

  outdir<-dirname(basefilename)
  basefile=file_path_sans_ext(basename(basefilename))

  cat('Analyzing',type,'\n')

  for(m in 1:length(algmodes)) {
    cat(algmodes[m],'\n')
    if( type == 'arel' || type == 'edge' || type == 'conf' ) {
      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-result.Rdata',sep=''),full.names=TRUE)
      if( length(resultfiles) >= 1 ) {
	load(resultfiles[1])
	positives<-matrix(0,result$p,result$p)
	negatives<-matrix(0,result$p,result$p)
        avg<-matrix(0,result$p,result$p)
	for(i in 1:length(resultfiles)) {
	  load(resultfiles[i])
	  if( type == 'arel' ) {
	    positives<-positives+(result$arel>0)
	    negatives<-negatives+(result$arel<0)
            avg<-avg+result$arel
	  } else if( type == 'edge' ) {
	    positives<-positives+(result$edge>0)
	    negatives<-negatives+(result$edge<0)
            avg<-avg+result$edge
	  } else if( type == 'conf' ) {
	    positives<-positives+(result$conf>0)
	    negatives<-negatives+(result$conf<0)
            avg<-avg+result$conf
          }
	}
	diff<-positives-negatives
	diff<-diff-diag(diag(diff))
	colnames(diff)<-result$labels
	cat(diff,'\n')
        avg<-avg / length(resultfiles)
	colnames(diff)<-result$labels
        cat(avg,'\n')
#	write.csv(diff,file=paste(basefilename,'-',algmodes[m],'-bs-',type,'.csv',sep=''),row.names=FALSE)
	write.csv(avg,file=paste(basefilename,'-',algmodes[m],'-bs-',type,'.csv',sep=''),row.names=FALSE)
      }
    }
    else if( type == 'runtime' ) {
      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*.runtime',sep=''),full.names=TRUE)
      total_runtime<-0
      if( length(resultfiles) >= 1 ) {
	for(i in 1:length(resultfiles)) {
	  timef<-file(resultfiles[i],'r')
	  runtime<-readLines(timef,1)
	  close(timef)
	  total_runtime <- total_runtime + as.numeric(runtime)
	}
	cat(file=paste(basefilename,'-',algmodes[m],'-bs.runtime',sep=''),total_runtime,'\n')
      }
    }
    else if( type=='cleanup' ) {
      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-result.Rdata',sep=''),full.names=TRUE)
      for(i in 1:length(resultfiles))
	unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*.runtime',sep=''),full.names=TRUE)
      for(i in 1:length(resultfiles))
	unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-arel.csv',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-arel.dot',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-edge.csv',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-edge.dot',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-conf.csv',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-pag.dot',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])

      resultfiles<-list.files(path=outdir,pattern=paste(basefile,'-',algmodes[m],'-[0123456789]*-mag.dot',sep=''),full.names=TRUE)
      if( length(resultfiles) > 0 )
	for(i in 1:length(resultfiles))
	  unlink(resultfiles[i])
    }
  }
}


args <- commandArgs(trailingOnly = TRUE)
basefilename <- args[1]
cleanup <- FALSE
if( !is.na(args[2]) )
  cleanup <- as.numeric(args[2])
algmodes_noconfedge<-c('icp-sc','icp-mc','cif-mc','cif-mccon','cif-mcsct','cif-mcsctcon','cif-sc','cif-sccon')
algmodes_confedge<-c('lcd-mc','lcd-mccon','lcd-mcsct','lcd-mcsctcon','lcd-sc','lcd-sccon','fci-obs','fci-pooled','fci-meta','fci-jci123','fci-jci123r','fci-jci1','fci-jci0')
analyze(basefilename,algmodes=c(algmodes_noconfedge,algmodes_confedge),type='arel')
analyze(basefilename,algmodes=algmodes_confedge,type='edge')
analyze(basefilename,algmodes=algmodes_confedge,type='conf')
analyze(basefilename,algmodes=c(algmodes_noconfedge,algmodes_confedge),type='runtime')
if( cleanup )
  analyze(basefilename,algmodes=c(algmodes_noconfedge,algmodes_confedge),type='cleanup')
