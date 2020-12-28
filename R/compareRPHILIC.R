#'Performs peak picking based on centWave algorithm and compares peak list in two chromatographies
#'
#'In the first step it performs peak picking on both HILIC and RP files according to centWave approach (See References for more details).
#'Afterwards, peak tables are splited on lists objects so that each element of the lists consists of chromatographic peaks clustered according to mass accuracy (threshold argument mzthr).
#'If doublecheck is set to TRUE; For each element of the list, mean mz is passed again to centWave as Region Of Interest (ROI) and in case extra peaks are detected they are added to the lists
#'After this Verifycation of peak picking, two columns entitled "numlist" and "chromatogr" are added. 
#'Numlist indicates the element number of the list (mz slice number) and "chromatogr" indicates chromatography which is either RP or HILIC.
#'Then, elements of both lists are compared and thus common peaks are detected and removed from the initial lists. Now initial lists are unique peaks in the chromatographies.
#'Again if doublecheck is set to TRUE; Unique peaks are passed as ROIs in centWave and if a pair peak is detected in the other chromatography,
#'the pair is moved to common peaks and is deleted from unique lists.
#'@param HILIC path where HILIC mzXML file is stored
#'@param RP path where RP mzXML file is stored
#'@param mzthr absolute mass accuracy error
#'@param doubleckeck TRUE or FALSE value as described in Description section
#'@param ... parameters of centWave peak picking algorithm (ppm, peakwidth, sn, prefilter, ...)
#'
#'@return Returns a list of 3 elements; common peaks, peaks unique in HILIC, peaks unique in RP
#'@references 
#'\url{http://msbi.ipb-halle.de/msbi/centwave} and DOI: 10.1186/1471-2105-9-504 
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 
compareRPHILIC<-function(HILIC=paste(getwd(),"Subtracted_HILIC.mzXML",sep="/"),
                  RP=paste(getwd(),"Subtracted_RP.mzXML",sep="/"),
                  mzthr=0.01, doubleckeck=TRUE, ...){ 

  #1. Peak Picking HILIC
  sample_HILIC_xset <- xcmsSet(HILIC, method = 'centWave',...)
  
  #2.Convert HILIC peak table matrix to a list object with elements having the same mass (Based on mass accuracy)
  peaktable_HILIC<-as.data.frame(sample_HILIC_xset@peaks)
  
  peaktable_HILIC<- peaktable_HILIC[order(peaktable_HILIC$mz),]
  listpeaktable_HILIC<-list()
  i<-1
  k<-1
  
  while(sum(peaktable_HILIC$mz)!=0){
    if(peaktable_HILIC$mz[i]!=0){
    if(min(abs(peaktable_HILIC$mz[i]-peaktable_HILIC$mz[-i]))<mzthr){
      listpeaktable_HILIC[[k]]<-rbind(peaktable_HILIC[abs(peaktable_HILIC$mz[i]-peaktable_HILIC$mz)<mzthr,])
    } else {
      listpeaktable_HILIC[[k]]<-peaktable_HILIC[i,]
    }
      peaktable_HILIC[abs(peaktable_HILIC$mz[i]-peaktable_HILIC$mz)<mzthr,]<-0
    k<-k+1
    }
    i<-i+1
  }
  ##Remove any duplicate peaks
  #i<-1
  #for(i in 1:length(listpeaktable_HILIC)){
  #  if(any(duplicated(listpeaktable_HILIC[[i]]$maxo))) listpeaktable_HILIC[[i]]<-listpeaktable_HILIC[[i]][duplicated(listpeaktable_HILIC[[i]]),]
  #}
  ##
  
  #3.Passes as ROI masses of each element of the list for re-peak picking 
  if(doubleckeck==T){
    progress<-txtProgressBar(min=1, max=length(listpeaktable_HILIC), style=3)
    i<-1
    peakHILIC<-c(0)
  for(i in 1:length(listpeaktable_HILIC)){ 
    j<-1
    capture.output(caseoferror<-try(
    peakHILIC<-xcmsSet(HILIC,method = 'centWave', ...,
                                                    ROI.list=list(data.frame(scmax=10000000,scmin=1,
                                                                             mzmin=c(mean(listpeaktable_HILIC[[i]]$mz)-mzthr),
                                                                             mzmax=c(mean(listpeaktable_HILIC[[i]]$mz)+mzthr))))
                                  ,silent=T))
  
  if(length(peakHILIC@peaks[,1])!=0){
  for(j in 1:length(peakHILIC@peaks[,1])){  
  if(!any(as.data.frame(peakHILIC@peaks)$maxo[j]==listpeaktable_HILIC[[i]]$maxo)){
    #print(i)
    #cat("In the slice","\n") 
    #print(as.data.frame(t(do.call(rbind.data.frame,listpeaktable_HILIC[[i]])),row.names = "1"))
    #cat("the following peak was added","\n")
    #print(data.frame(t(as.data.frame(peakHILIC@peaks[j,])),row.names = "1"))
    listpeaktable_HILIC[[i]][c(length(listpeaktable_HILIC[[i]]$mz)+1),]<-peakHILIC@peaks[j,]
  }
                                         }
  }
    setTxtProgressBar(progress, i)
  }
    cat("\n","Verifying peak picking is done...","\n")
}
##############################
  
  #1. Peak Picking RP
  sample_RP_xset <- xcmsSet(RP,method = 'centWave', ...)

  
  #2.Convert HILIC peak table matrix to a list object with elements having the same mass (Based on mass accuracy)
  peaktable_RP<-as.data.frame(sample_RP_xset@peaks)
  peaktable_RP<- peaktable_RP[order(peaktable_RP$mz),]
  listpeaktable_RP<-list()
  i<-1
  k<-1
  
  while(sum(peaktable_RP$mz)!=0){
    if(peaktable_RP$mz[i]!=0){
      if(min(abs(peaktable_RP$mz[i]-peaktable_RP$mz[-i]))<mzthr){
        listpeaktable_RP[[k]]<-rbind(peaktable_RP[abs(peaktable_RP$mz[i]-peaktable_RP$mz)<mzthr,])
      } else {
        listpeaktable_RP[[k]]<-peaktable_RP[i,]
      }
      peaktable_RP[abs(peaktable_RP$mz[i]-peaktable_RP$mz)<mzthr,]<-0
      k<-k+1
    }
    i<-i+1
  }
  
  ##Remove any duplicate peaks
  #i<-1
  #for(i in 1:length(listpeaktable_RP)){
  #  if(any(duplicated(listpeaktable_RP[[i]]$maxo))) listpeaktable_RP[[i]]<-listpeaktable_RP[[i]][duplicated(listpeaktable_RP[[i]]),]
  #}
  ##
  
  #3.Passes as ROI masses of each element of the list for re-peak picking
  if(doubleckeck==T){
  progress<-txtProgressBar(min=1, max=length(listpeaktable_RP), style=3)
  i<-1
  for(i in 1:length(listpeaktable_RP)){ 
    j<-1
    capture.output(caseoferror<-try(
      peakRP<-xcmsSet(RP,method = 'centWave', ...,
                         ROI.list=list(data.frame(scmax=10000000,scmin=1,
                                                  mzmin=c(mean(listpeaktable_RP[[i]]$mz)-mzthr),
                                                  mzmax=c(mean(listpeaktable_RP[[i]]$mz)+mzthr))))
      ,silent=T))
    if(length(peakRP@peaks[,1])!=0){
      for(j in 1:length(peakRP@peaks[,1])){  
        if(!any(as.data.frame(peakRP@peaks)$maxo[j]==listpeaktable_RP[[i]]$maxo)){
          #print(i)
          #cat("In the slice","\n") 
          #print(as.data.frame(t(do.call(rbind.data.frame,listpeaktable_RP[[i]])),row.names = "1"))
          #cat("the following peak was added","\n")
          #print(data.frame(t(as.data.frame(peakRP@peaks[j,])),row.names = "1"))
          listpeaktable_RP[[i]][c(length(listpeaktable_RP[[i]]$mz)+1),]<-peakRP@peaks[j,]
        }
      }
    }
    setTxtProgressBar(progress, i)
  }
  cat("\n","Verifying peak picking on the other sample is done...","\n")
  }
  
  #4.Add numlist and chromatogr in the peak lists
  i<-1
  for(i in 1:length(listpeaktable_RP)){
  listpeaktable_RP[[i]]$numlist<-i
  listpeaktable_RP[[i]]$chromatogr<-"RP"
  }
  #####
  i<-1
  for(i in 1:length(listpeaktable_HILIC)){
    listpeaktable_HILIC[[i]]$numlist<-i
    listpeaktable_HILIC[[i]]$chromatogr<-"HILIC"
  }

  #5.Get common slices  
  peaktable_RP2<-do.call(rbind.data.frame,listpeaktable_RP)
  peaktable_HILIC2<-do.call(rbind.data.frame,listpeaktable_HILIC)
  pairs<-data.frame(HILIC=0,RP=0)
  i<-1
  k<-1
  for(i in 1:length(peaktable_HILIC2[,1])){
  if(min(abs(peaktable_HILIC2[i,1]-peaktable_RP2[,1]))<mzthr){
    temp<-c(peaktable_HILIC2$numlist[i],peaktable_RP2$numlist[which.min(abs(peaktable_HILIC2[i,1]-peaktable_RP2[,1]))])
    if(k==1 || pairs[k-1,]!=temp){
    pairs[k,]<-temp
    k<-k+1
    }
      }
  }
  #####
  peaktable_RP2<-do.call(rbind.data.frame,listpeaktable_RP)
  peaktable_HILIC2<-do.call(rbind.data.frame,listpeaktable_HILIC)
  pairs2<-data.frame(RP=0,HILIC=0)
  i<-1
  k<-1
  for(i in 1:length(peaktable_RP2[,1])){
    if(min(abs(peaktable_RP2[i,1]-peaktable_HILIC2[,1]))<mzthr){
     # print(i)
      temp<-c(peaktable_RP2$numlist[i],peaktable_HILIC2$numlist[which.min(abs(peaktable_RP2[i,1]-peaktable_HILIC2[,1]))])
      if(k==1 || pairs2[k-1,]!=temp){
        pairs2[k,]<-temp
        k<-k+1
      }
    }
  }
  #####
  pairs$concatenated<-as.numeric(paste(pairs[,1],pairs[,2], sep=""))
  pairs2$concatenated<-as.numeric(paste(pairs2[,2],pairs2[,1], sep=""))
  i<-1
  for(i in 1:length(pairs2$concatenated)){
    if(!any(pairs2$concatenated[i]==pairs$concatenated)){
       pairs[c(length(pairs[,1])+1),]<-pairs[i,]
      }
  }
  

  i<-1
  common<-list()
  for(i in 1:length(pairs[,1])){
  common[[i]]<-rbind(listpeaktable_HILIC[[pairs[i,1]]],listpeaktable_RP[[pairs[i,2]]])
  }

  uniqueHILIC<-listpeaktable_HILIC[-pairs[,1]]
  uniqueRP<-listpeaktable_RP[-pairs[,2]]


  i<-1
  for(i in 1:length(uniqueHILIC)){
  uniqueHILIC[[i]]$numlist<-i
  uniqueHILIC[[i]]$chromatogr<-"HILIC"
  }
  i<-1
  for(i in 1:length(uniqueRP)){
  uniqueRP[[i]]$numlist<-i
  uniqueRP[[i]]$chromatogr<-"RP"
  }

  
  #6.For unique peaks in HILIC, we pass mz as ROI for peak picking in RP. If a peak exists in RP then the pair is moved to common list. 
  i<-1
  progress<-txtProgressBar(min=1, max=length(uniqueHILIC), style=3)
  k<-1
  lc<-c(length(common))
  remove_lines<-c(rep(TRUE,length(uniqueHILIC)))
  for(i in 1:length(uniqueHILIC)){
    capture.output(caseoferror<-try(pairpeak<-xcmsSet(RP,method = 'centWave', ...,
    ROI.list=list(data.frame(scmax=10000000,scmin=1,mzmin=c(uniqueHILIC[[i]]$mz-mzthr),mzmax=c(uniqueHILIC[[i]]$mz+mzthr)))),silent=T))
    if(!is.character(caseoferror)){
      if(length(pairpeak@peaks[,1])!=0){
       # print(i)
        remove_lines[i]<-FALSE
        pairpeak<-as.data.frame(pairpeak@peaks)
        pairpeak$numlist<-"NO"
        pairpeak$chromatogr<-"RP"
        backtocommonlist<-rbind(pairpeak,uniqueHILIC[[i]])
        common[[c(lc+k)]]<-backtocommonlist
        k<-k+1
      }
    }
    setTxtProgressBar(progress, i)
  }
  cat("\n","Crossverifying of unique peaks is done...","\n")
  uniqueHILIC<-uniqueHILIC[remove_lines]

  
  
  
  
  #7.For unique peaks in RP, we pass mz as ROI for peak picking in HILIC. If a peak exists in HILIC then the pair is moved to common list. 
  i<-1
  progress<-txtProgressBar(min=1, max=length(uniqueRP), style=3)
  k<-1
  lc<-c(length(common))
  remove_lines<-c(rep(TRUE,length(uniqueRP)))
  for(i in 1:length(uniqueRP)){
    capture.output(caseoferror<-try(pairpeak<-xcmsSet(HILIC,method = 'centWave', ...,
                                                      ROI.list=list(data.frame(scmax=10000000,scmin=1,mzmin=c(mean(uniqueRP[[i]]$mz)-10*mzthr),mzmax=c(mean(uniqueRP[[i]]$mz)+10*mzthr)))),silent=T))
    if(!is.character(caseoferror)){
      if(length(pairpeak@peaks[,1])!=0){
       # print(i)
        remove_lines[i]<-FALSE
        pairpeak<-as.data.frame(pairpeak@peaks)
        pairpeak$numlist<-"NO"
        pairpeak$chromatogr<-"HILIC"
        backtocommonlist<-rbind(pairpeak,uniqueRP[[i]])
        common[[c(lc+k)]]<-backtocommonlist
        k<-k+1
      }
    }
    setTxtProgressBar(progress, i)
  }
  cat("\n","Crossverifying of unique peaks in the second sample is done...","\n")
  uniqueRP<-uniqueRP[remove_lines]

  
  output<-list()
  output[[1]]<-common
  output[[2]]<-uniqueHILIC
  output[[3]]<-uniqueRP
  
  
  return(output)
}
