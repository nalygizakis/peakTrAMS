#'Extracts and saves MS/MS spectra from a directory with mzXML files 
#'
#'It is possible (especially in negative ionization) some MS/MS spectrum to be empty (without any spectral peaks).
#'This function removes the reduntant full scans and returns a mzXML object without the empty full scans.
#'@param sample mzXML list produced by read.mzXML function
#'@return sample mzXML list without empty full scans
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'
#'@export 

extractMS2<-function(mzXML){
files<-list.files("E:/3.Treatment_chain/151019_RP pos/mzXML",full.names = T)
files_short<-list.files("E:/3.Treatment_chain/151019_RP pos/mzXML",full.names = F)
parent_directory<-getwd()

k<-1
for(k in 1:length(files_short)){
setwd(parent_directory)
filemzxml<-read.mzXML(files[k])
dir.create(paste(files_short[k]))
setwd(paste(files_short[k]))
info_filemzxml<-getinfo(filemzxml)
info_filemzxml<-info_filemzxml[info_filemzxml$mslevel==2,]
info_filemzxml<-info_filemzxml[info_filemzxml$numofpeaks>0,]
info_filemzxml<-info_filemzxml[info_filemzxml$timeofscan>15.6,]

i<-1
for(i in 1:length(info_filemzxml[,1])){
peaklist<-cbind(round(filemzxml$scan[[info_filemzxml[i,1]]][[1]],4),filemzxml$scan[[info_filemzxml[i,1]]][[2]])
namepeaklist<-paste(round(as.numeric(strsplit(strsplit(filemzxml$scan[[info_filemzxml[i,1]]]$precursorMz," </precursorMz")[[1]][1],"\">   ")[[1]][2]),4),
                    round(as.numeric(strsplit(strsplit(filemzxml$scan[[info_filemzxml[i,1]]]$scanAttr,"PT")[[1]][2],"S\"")[[1]][1])/60,2),sep="_")
namepeaklist<-paste(namepeaklist,"RP","POS",sep="_")
peaklist<-peaklist[peaklist[,1]<c(round(as.numeric(strsplit(strsplit(filemzxml$scan[[info_filemzxml[i,1]]]$precursorMz," </precursorMz")[[1]][1],"\">   ")[[1]][2]),4)+2.5),]

write.table(peaklist,file=paste(namepeaklist,".txt"),row.names=F,sep=" ",col.names = FALSE)
}
print(k)
}
return(NULL)
}