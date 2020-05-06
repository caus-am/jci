# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

dg2graphviz <- function(filename,p,contextVars,labels,dg) {
  sink(filename)
  cat('digraph G {\n');
  cat('subgraph cluster_system {\n');
  cat('label="system";\n')
  cat('style="invis";\n')
  if( p > 0 ) {
    for( i in 1:p )
      if( !(i %in% contextVars) )
        cat(i,'[label="',labels[i],'"];\n',sep='')
    for( i in 1:p )
      if( !(i %in% contextVars) ) {
        for( j in 1:p )
          if( !(j %in% contextVars) && i!=j )
            if( dg[i,j] ) {
              cat(i,'->',j,';\n',sep='')
            }
      }
  }
  cat('}\n')
  cat('subgraph context {\n');
  cat('label="context";\n')
  if( p > 0 ) {
    for( i in 1:p )
      if( i %in% contextVars )
        cat(i,'[label="',labels[i],'",shape=rectangle];\n',sep='')
    for( i in 1:p )
      if( i %in% contextVars ) {
        for( j in 1:p )
          if( !(j %in% contextVars) ) {
            if( dg[i,j] ) 
              cat(i,'->',j,';\n',sep='')
          }
      }
  }
  cat('}\n')
  cat('}\n')
  sink()
}
