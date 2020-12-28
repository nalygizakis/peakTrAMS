#'Removes all MS2 scan events from a mzXML list
#'
#'Takes as input an object which was created by read.mzXML function 
#'and returns an object without the MS/MS full scans.
#'@param sample mzXML list object
#'@return sample mzXML list object without MS/MS spectra
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 
removeMS2<-function(sample){
  
  i<-1; removescanevents<-c();
  for(i in 1:length(sample$scan)) removescanevents[i]<-sample$scan[i][[1]]$msLevel
  
  removescanevents2<-c(); i<-1;
  for(i in 1:length(removescanevents))  if((removescanevents[i]-1)) removescanevents2[i]<-i
  

  new_sample<-list()
  new_sample<-sample[1:4]
  new_sample$scan<-sample[[5]][-removescanevents2[!is.na(removescanevents2)]]

  attr(new_sample, "class") = "mzXML"
  
  i<-1
  for(i in 1:length(new_sample$scan)) new_sample[[5]][[i]]$num<-i
  
  return(new_sample)
}