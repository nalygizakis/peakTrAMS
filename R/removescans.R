#'Removes selected full scans from a mzXML list object
#'
#'Takes as input an object which was created by read.mzXML function 
#'and returns an object without selected scans passed in scansORtime argument.
#'In case time is set to TRUE then scanORtime should be a vector of two elements containing
#'retention time in minutes. The selected full scans with retention time within this interval will
#'be deleted from the mzXML list.
#'@param mzXML file produced from read.mzXML function
#'@param scansORtime Selected scans (or scans with retention time if time=TRUE) to be removed
#'@param time Logical value TRUE or FALSE
#'@return a mzXML list without selected full scans
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export
removescans<-function(mzXML=blank_HILIC,scansORtime=c(17,25),time=TRUE){

  if(scansORtime[2]=="end" & time==FALSE) scansORtime[2]<-max(getinfo(mzXML)$scan)
  if(scansORtime[2]=="end" & time==TRUE) scansORtime[2]<-max(getinfo(mzXML)$timeofscan)/60
  if(scansORtime[1]=="beginning" & time==FALSE) scansORtime[1]<-min(getinfo(mzXML)$scan)
  if(scansORtime[1]=="beginning" & time==TRUE) scansORtime[1]<-min(getinfo(mzXML)$timeofscan)/60
  scansORtime<-as.numeric(scansORtime)
  scansORtime2<-scansORtime
  
  info<-getinfo(mzXML)
  
  if(time==TRUE){
  k<-which.min(abs(info$timeofscan-scansORtime[2]*60))
  if(info$mslevel[which.min(abs(info$timeofscan-scansORtime[2]*60))]!=1 & k!=length(info[,1])){
  while(info$mslevel[k]!=1) {
    k <- k+1
    scansORtime[2]<-info$timeofscan[k]/60
  }
  k<-k-1
  cat("Ending point was set at", paste(round(c(info$timeofscan[k]/60),4)), "because given ending retention time", scansORtime2[2] ,"corresponds to scan at MS2 level","\n")
  }
  } else {
    #Prevent removal of MSn scans
    scansORtime<-scansORtime[info$mslevel[scansORtime]==1]
  }
  
  if(time==TRUE){
  u<-which.min(abs(info$timeofscan-scansORtime[1]*60))
  if(info$mslevel[which.min(abs(info$timeofscan-scansORtime[1]*60))]!=1){
    info<-getinfo(mzXML)
    while(info$mslevel[u]!=1) {
      u <- u-1 
      scansORtime[1]<-info$timeofscan[u]/60
    }
    u<-u-1
    cat("Beginning point was set at", paste(round(c(info$timeofscan[u]/60),4)), "because given ending retention time", scansORtime2[1] ,"corresponds to scan at MS2 level","\n")
  }
  } else {
    #Prevent removal of MSn scans
    scansORtime<-scansORtime[info$mslevel[scansORtime]==1]
  }
  
  if(time==TRUE){
  if(!is.na(info$mslevel[which.min(abs(info$timeofscan-scansORtime[2]*60))+1]!=1)){
  if(info$mslevel[which.min(abs(info$timeofscan-scansORtime[2]*60))]==1 & info$mslevel[which.min(abs(info$timeofscan-scansORtime[2]*60))+1]!=1) k<-which.min(abs(info$timeofscan-scansORtime[2]*60))-1
  }
  if(!is.na(info$mslevel[which.min(abs(info$timeofscan-scansORtime[1]*60))+1]!=1)){
  if(info$mslevel[which.min(abs(info$timeofscan-scansORtime[1]*60))]==1 & info$mslevel[which.min(abs(info$timeofscan-scansORtime[1]*60))+1]!=1) u<-which.min(abs(info$timeofscan-scansORtime[1]*60))-1
  }
  }
  #Output of the procedure is scansORtime object
  
  if(time==TRUE){
    info<-getinfo(mzXML)[u:k,]
    stayordelete<-c(rep(TRUE,length(mzXML$scan)))
    stayordelete[info$scan]<-FALSE
    } else {
      stayordelete<-c(rep(TRUE,length(mzXML$scan)))
      stayordelete[scansORtime]<-FALSE
    }
  
  new_sample<-list()
  new_sample<-mzXML[1:4]
  new_sample$scan<-mzXML$scan[c(stayordelete)]
  attr(new_sample, "class") = "mzXML"
  
  
  i<-1
  for(i in 1:length(new_sample$scan)) new_sample[[5]][[i]]$num<-i
  
  return(new_sample)
  }