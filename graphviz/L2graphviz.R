# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

L2graphviz <- function(filename,p,contextVars,labels,L,mode='jci123') {
   # The output L is a model as used in the HEJ2014 code
   # it has fields G, Ge, Gs, Gcircles
   #   G are the directed edges
   #   Ge the bidirected edges
   #   Gs undirected edges
   #   Gcircles are unused
   # e.g. G[to,from]=1 iff from --> to (i.e. the opposite convention as the pag adjacency matrix in pcalg)
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
          if( !(j %in% contextVars) ) {
            if( L$G[i,j] ) 
              cat(j,'->',i,';\n',sep='')
            if( L$G[j,i] )
              cat(i,'->',j,';\n',sep='')
            if( L$Ge[i,j] )
              cat(i,'->',j,'[dir="both"];\n',sep='')
            if( L$Gs[i,j] )
              cat(i,'->',j,'[arrowhead="none"];\n',sep='')
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
              if( L$G[i,j] ) 
                cat(j,'->',i,';\n',sep='')
              if( L$G[j,i] ) 
                cat(i,'->',j,';\n',sep='')
              if( L$Ge[i,j] )
                cat(i,'->',j,'[dir="both"];\n',sep='')
              if( L$Gs[i,j] )
                cat(i,'->',j,'[arrowhead="none"];\n',sep='')
            }
        }
    } else {
      if( p >= 2 ) {
        for( i in 1:(p-1) )
          for( j in (i+1):p )
            if( (i %in% contextVars) || (j %in% contextVars) ) {
              if( L$G[i,j] ) 
                cat(j,'->',i,';\n',sep='')
              if( L$G[j,i] ) 
                cat(i,'->',j,';\n',sep='')
              if( L$Ge[i,j] )
                cat(i,'->',j,'[dir="both"];\n',sep='')
              if( L$Gs[i,j] )
                cat(i,'->',j,'[arrowhead="none"];\n',sep='')
            }
      }
    }
  }
#  cat('}\n')
  cat('}\n')
  sink()
}
