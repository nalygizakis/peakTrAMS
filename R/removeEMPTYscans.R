#'Removes empty scan events from a mzXML file
#'
#'It is possible (especially in negative ionization) some MS/MS spectrum to be empty (without any spectral peaks).
#'This function removes the reduntant full scans and returns a mzXML object without the empty full scans.
#'@param sample mzXML list produced by read.mzXML function
#'@return new_sample mzXML list without empty full scans
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 

removeEMPTYscans<-function(sample){
    info<-getinfo(sample)
    remove<-info[info$numofpeaks==0,]

    new_sample<-list()
    new_sample<-sample[1:4]
    new_sample$scan<-sample[[5]][-remove$scan]
    attr(new_sample, "class") = "mzXML"
    
    i<-1
    for(i in 1:length(new_sample$scan)) new_sample[[5]][[i]]$num<-i
    
    return(sample)
    }
 