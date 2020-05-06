# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

pag2graphviz <- function(filename,p,contextVars,labels,pag,mode='jci123') {
  sink(filename)
  cat('digraph G {\n');
#  cat('subgraph cluster_system {\n');
#  cat('label="system";\n')
#  cat('style="invis";\n')
  arrowtypes=c("odot","normal","none")
  if( p > 0 ) {
    for( i in 1:p )
      if( !(i %in% contextVars) )
        cat(i,'[label="',labels[i],'"];\n',sep='')
  }
  if( p >= 2 ) {
    for( i in 1:(p-1) )
      if( !(i %in% contextVars) ) {
        for( j in (i+1):p )
          if( !(j %in% contextVars) )
            if( pag[i,j] ) {
              cat(i,'->',j,'[dir="both",arrowhead="',arrowtypes[pag[i,j]],'",arrowtail="',arrowtypes[pag[j,i]],'"];\n',sep='')
            }
      }
  }
#  cat('}\n')
#  cat('subgraph context {\n');
#  cat('label="context";\n')
  if( p > 0 ) {
    for( i in 1:p )
      if( i %in% contextVars )
        cat(i,'[label="',labels[i],'",shape=rectangle];\n',sep='')
    if( mode == 'jci123' ) {
      cat('R[shape=rectangle,style="dashed"];\n',sep='')
      for( i in 1:p )
	if( i %in% contextVars ) {
	  cat('R->',i,'[style="dashed"];\n',sep='')
	  for( j in 1:p )
	    if( !(j %in% contextVars) ) {
	      if( pag[i,j] ) 
		cat(i,'->',j,'[dir="both",arrowhead="',arrowtypes[pag[i,j]],'",arrowtail="',arrowtypes[pag[j,i]],'"];\n',sep='')
	    }
	}
    } else {
      if( p >= 2 ) {
        for( i in 1:(p-1) )
          for( j in (i+1):p )
            if( (i %in% contextVars) || (j %in% contextVars) )
              if( pag[i,j] ) {
                cat(i,'->',j,'[dir="both",arrowhead="',arrowtypes[pag[i,j]],'",arrowtail="',arrowtypes[pag[j,i]],'"];\n',sep='')
              }
      }
    }
  }
#  cat('}\n')
  cat('}\n')
  sink()
}
