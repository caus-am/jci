# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

datafile<-'../../../../../aci/experiments/sachs/data/sachs.csv'
data<-read.csv(file=datafile,header=TRUE,sep="\t")
augdatafile<-'sachs-jci-noICAM.csv'

expdesign<-matrix(0,8,6)
 expdesign[1 ,] <- c(      0, 0, 0, 0, 0, 0)
#expdesign[2 ,] <- c(   1, 0, 0, 0, 0, 0, 0)
 expdesign[2 ,] <- c(      1, 0, 0, 0, 0, 0)
 expdesign[3 ,] <- c(      0, 1, 0, 0, 0, 0)
 expdesign[4 ,] <- c(      0, 0, 1, 0, 0, 0)
 expdesign[5 ,] <- c(      0, 0, 0, 1, 0, 0)
 expdesign[6 ,] <- c(      0, 0, 0, 0, 1, 0)
 expdesign[7 ,] <- c(      0, 0, 0, 0, 0, 1)
 expdesign[8 ,] <- c(      0, 0, 0, 0, 0, 2)
#expdesign[10,] <- c(   1, 1, 0, 0, 0, 0, 0)
#expdesign[11,] <- c(   1, 0, 1, 0, 0, 0, 0)
#expdesign[12,] <- c(   1, 0, 0, 1, 0, 0, 0)
#expdesign[13,] <- c(   1, 0, 0, 0, 1, 0, 0)
#expdesign[14,] <- c(   1, 0, 0, 0, 0, 1, 0)

newR<-c(1,NA,2,3,4,5,6,7,8,NA,NA,NA,NA,NA)
Cindex<-dim(data)[2]
augdata<-matrix(0,dim(data)[1],Cindex+dim(expdesign)[2])
row<-0
for( orgrow in 1:(dim(data)[1]) ) {
  R<-data[orgrow,Cindex]
  if( R %in% which(!is.na(newR)) ) {
    row<-row+1
    augdata[row,]<-c(as.numeric(data[orgrow,1:(Cindex-1)]),R,expdesign[newR[R],])
  }
}

augdata<-data.frame(augdata[1:row,])
colnames(augdata)<-c(colnames(data),"AKT inh","G0076","Psitectorigenin","U0126","LY294002","PMA/beta2CAMP + noAlphaCD3/28")
augdata[,1:11]<-log(augdata[,1:11])
write.csv(augdata,file=augdatafile,row.names=FALSE)
