#'Gets retention time and number of peaks of full scans of a mzXML list
#'
#'Takes in a raw sample and returns a data frame with retention time of each full scan
#'@param sample mzXML list created from read.mzXML function
#'@return A data frame with number of scan, with retention time of each full scan, mslevel, number of spectral peaks and in case
#'of MS/MS full scan precursor mass and precurson intensity.
#'
#'@examples
#'sample_mzXML<-read.mzXML(list.files(paste(find.package(package="peakTrams"),"data",sep="/"),pattern = ".mzXML", full.names = TRUE))
#'getrt(sample=sample_mzXML)
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'
#'@export
getinfo<-function(sample){
  numscan<-sample$scan[[1]]$num

  info<-data.frame(scan=numscan:length(sample$scan),timeofscan=0)
  for(numscan in 1:c(length(sample$scan)-numscan+1)){ 
    if(length(strsplit(try(sample$scan[[numscan]][[6]], silent=T),"Error")[[1]])!=2){
    info$timeofscan[numscan]<-sample$scan[[numscan]][[6]]
    info[numscan,2]<-as.numeric(strsplit(strsplit(info[numscan,2],split="S")[[1]][1],split="PT")[[1]][2])
  #  info$basePeakMz[numscan]<-sprintf("%.5f",sample$scan[[i]]$mass[which.max(sample$scan[[numscan]]$mass)])
  #  info$basePeakIntensity[numscan]<-as.numeric(sprintf("%.0f",max(sample$scan[[numscan]]$peaks)))
    }
  }
  info$timeofscan<-as.numeric(info$timeofscan)
  
  numscan<-sample$scan[[1]]$num
  for(numscan in 1:c(length(sample$scan)-numscan+1)){ 
    if(length(strsplit(try(sample$scan[[numscan]][[6]], silent=T),"Error")[[1]])!=2){
    info$mslevel[numscan]<-(sample$scan[[numscan]][[5]])
    info$numofpeaks[numscan]<-length(sample$scan[[numscan]][[1]])
    info$CE[numscan]<-strsplit(strsplit(sample$scan[[numscan]]$scanAttr, "collisionEnergy=[\"]")[[1]][2], "[\"]")[[1]][1]
    }
  }
  info$precursor<-NA
  info$precursorIntensity<-NA
  
  i<-1
  for(i in 1:length(info[,1])){
    if(length(strsplit(try(sample$scan[[numscan]][[6]], silent=T),"Error")[[1]])!=2){
  if(info$mslevel[i]!=1){
  info$precursor[i]<-as.numeric(strsplit(strsplit(sample$scan[[i]]$precursorMz," </precursorMz>\n")[[1]][1],">   ")[[1]][2])
  info$precursorIntensity[i]<-as.numeric(strsplit(sample$scan[[i]]$precursorMz,"[\"]")[[1]][2])
    }
    }
  }
  sprintf("Done")
  
  info<-info[info$timeofscan!=0,]
  return(info)
}