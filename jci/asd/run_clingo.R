# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

run_clingo <- function(basefilename,extrafiles,query='',clingo_cmd='clingo',clingo_pars='--quiet=2,1 -W no-atom-undefined',verbose=0) {
  suppressMessages(library("tools"))

  basedir=dirname(basefilename)
  basefile=file_path_sans_ext(basename(basefilename))
  baseext=file_ext(basefilename)

  if( query != '' ) {
    # write query to temporary file
    queryfile<-tempfile(pattern=basefile,tmpdir=basedir, fileext=".query")
    cat(file=queryfile,query,'\n')
  } else
    queryfile<-''

  outfile<-tempfile(pattern=basefile,tmpdir=basedir,fileext=".out")
  # run clingo
  cmd=paste(clingo_cmd,clingo_pars,basefilename)
  if( length(extrafiles) > 0 )
    for( i in 1:length(extrafiles) )
      cmd=paste(cmd,extrafiles[i])
  cmd=paste(cmd,queryfile,'>',outfile)
  if( verbose ) 
    cat('Running',cmd,'\n')
  status<-system(cmd)

  # read and remove output file
  outf<-file(outfile,'r')
  solution<-readLines(outf,-1)
  close(outf)
  unlink(outfile)

  if( verbose ) {
    cat('---\n')
    cat(solution,sep='\n')
    cat('---\n')
  }

  # parse output
  if( length(grep('^UNSATISFIABLE$',solution)) > 0 )
    val<-Inf
  else if( length(grep('^SATISFIABLE$',solution)) > 0 )
    val<-0
  else {
    result <- grep('^Optimization : ',solution)
    if( length(result) != 1 ) {
      cat('ERROR: ',cmd,'\n')
      stop(paste('clingo output for',queryfile,'not as expected',status))
    }
    # get third word from that line
    val<-as.numeric(strsplit(solution[result],' ')[[1]][3])
  }
  if( verbose )
    cat('Value: ',val,'\n')

  if( query != '' ) {
    # remove temporary query file
    unlink(queryfile)
  }
  
  return(val)
}
