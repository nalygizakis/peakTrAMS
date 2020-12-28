#'It searches MS/MS spectra for a mass within mass and retention time tolerance
#'
#'The function searches MS/MS full scans and looks for the selected mass within the given
#'tolerance mass. The MS/MS spectrum should be within a given retention time window, as this
#'is specified in rt and tolerancert arguments. In case many MS/MS spectrum exist in the file, the one with the most peak is selected. 
#'@param sample mzXML list produced by read.mzXML
#'@param mass selected mass
#'@param tolerancemass absolute mass accuracy set as tolerance value
#'@param rt retention time in minutes
#'@param tolerancert tolerance in retention time (given in minutes)
#'
#'@return Peak list data frame (mz and intensity) with MS/MS spectra. 
#'
#'@examples
#'sample_mzXML<-read.mzXML(list.files(paste(find.package(package="peakTrams"),"data",sep="/"),pattern = ".mzXML", full.names = TRUE))
#'searchMS2(sample=sample_mzXML, mass=100.0323, tolerancemass=0.05, rt=3.00, tolerancert=0.5)
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 
searchMS2<-function(sample, mass=100.0323, tolerancemass=0.05, rt=3.00, tolerancert=0.5){
  
  rtsample<-getinfo(sample)
  rtsample<-rtsample[rtsample$mslevel==2,]
  
  i<-1
  for(i in 1:length(rtsample[,1])){
  rtsample$precursor[i]<-as.numeric(strsplit(strsplit(sample[[5]][[rtsample[i,1]]]$precursorMz ,">")[[1]][2] ,"</precursorMz")[[1]][1])
  }
  
  
  if(min(abs(rtsample$precursor-mass))<=tolerancemass){
  #cat(length(rtsample[abs(rtsample$precursor-mass)<tolerancemass,1]),"scans were detected within given mass tolerance")  
    temp<-rtsample[abs(rtsample$precursor-mass)<tolerancemass,]
    temp<-temp[c(temp$timeofscan/60)<c(rt+tolerancert) & c(temp$timeofscan/60)>c(rt-tolerancert),]
    if(length(temp[,1])==0){
      cat("No MS/MS spectra found")
      peaklist<-NULL
    } else  {
      peaklist<-cbind(sample[[5]][[temp[which.max(temp$numofpeaks),1]]]$mass,sample[[5]][[temp[which.max(temp$numofpeaks),1]]]$peaks)
      colnames(peaklist)<-c("mz","int")
      peaklist<-peaklist[peaklist[,1]<mass+tolerancemass,]
          }
  

    return(peaklist)
  }
}
