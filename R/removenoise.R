#'Removes noisy spectral peaks below an intensity threshold in all full scans
#'
#'The function goes in all scans of a specific mslevel, reads the spectra and removes peaks below noisethreshold.
#'@param mzXML mzXML list produced by read.mzXML function
#'@param noisethreshold noise threshold
#'@param mslevel ms level, can be c(1) or c(2) or c(1,2)
#'
#'@return Returns a new mzXML list without spectral peaks below noisethreshold
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export
removenoise<-function(mzXML=read.mzXML(list.files(paste(find.package(package="peakTrams"),"data",sep="/"),pattern = ".mzXML", full.names = TRUE)), 
                      noisethreshold=100, mslevel=c(1)){
info_mzXML<-getinfo(mzXML)
info_mzXML<-info_mzXML[info_mzXML$mslevel==mslevel,]

i<-1
for(i in 1:length(info_mzXML[,1])){
which<-mzML$scan[[info_mzXML[i,1]]][[2]]>noisethreshold
mzML$scan[[info_mzXML[i,1]]][[2]]<-mzML$scan[[info_mzXML[i,1]]][[2]][mzML$scan[[info_mzXML[i,1]]][[2]]>noisethreshold]
mzML$scan[[info_mzXML[i,1]]][[1]]<-mzML$scan[[info_mzXML[i,1]]][[1]][which]



}
return(mzXML)
}
