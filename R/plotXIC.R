#'Performs eΧtracted Ιon Chromatogram (XIC)
#'
#'Performs extracted ion chromatogram using specific mass accuracy window and optionally retention time window. 
#'In case mzrange is a vector with one element, accuracy argument should be specified, otherwise mzrange should be a vector 
#'of two elements. 
#'Extra title can be a vector with two elements. The first element is added in the beginning of the title and the second at the end of
#'the title.
#'@param sample mzXML list object (e.g. read.mzXML("c:/myfolder/sample.mzXML")
#'@param mzrange mz
#'@param accuracy absolute mass accuracy
#'@param timerange time range 
#'@param ploteic If TRUE plots extracted ion chromatogram otherwise just returns the dataframe object without plotting
#'@param type what type of plot should be drawn. Possible types are
#'@param extra_title Extra title in the beginning of the main title in figure
#'@param ... other arguments passed to plot function (such as type for type of line,col for color of line etc)
#'
#'@return A dataframe with three columns mz,intensity and timeofscan
#'
#'@examples plotXIC(sample=paste(list.files(file.path(find.package("peakTrams"), "data"),pattern=c("Sample_RP",".mzXML"),full.names = T)), mzrange=c(286.1438), accuracy=0.005, type="l")
#'@examples plotXIC(sample=paste(list.files(file.path(find.package("peakTrams"), "data"),pattern=c("Sample_RP",".mzXML"),full.names = T)), mzrange=c(286.1438), accuracy=0.005, type="l", col="red", timerange=c(2,3.5))
#'@examples plotXIC(sample=paste(list.files(file.path(find.package("peakTrams"), "data"),pattern=c("Sample_RP",".mzXML"),full.names = T)), mzrange=c(286.1438), accuracy=0.005, type="l", col="red", timerange=c(2,3.5), extra_title=c("EIC of","Morphine [M+H]"))
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'
#'@export 
plotXIC<-function(sample=sub_RP,
                  mzrange=c(272.0275), accuracy=0.005, 
                  timerange=c(round(min(getinfo(sample)$timeofscan/60),1),round(max(getinfo(sample)$timeofscan/60),1)),
                  ploteic=TRUE, extra_title="EIC of",...){
  
  timepoints<-getinfo(sample)
  timepoints$element<-1:length(timepoints[,1])
  timerange_old<-c(round(min(timepoints$timeofscan/60),1),round(max(timepoints$timeofscan/60),1))
  if(length(extra_title)==1) extra_title[2]<-""
  if(length(mzrange)==1) mzrange=c(mzrange-accuracy,mzrange+accuracy)
  timepoints<-timepoints[timepoints$mslevel==1,]
  timepoints$timeofscan<-c(timepoints$timeofscan/60)
  timepoints<-timepoints[timepoints$timeofscan<timerange[2],]
  timepoints<-timepoints[timepoints$timeofscan>timerange[1],]

  
  i<-1
  eic<-data.frame(mz=c(rep(0,times=length(timepoints[,1]))),int=c(rep(0,times=length(timepoints[,1]))),timeofscan=c(rep(0,times=length(timepoints[,1]))))
  for(i in 1:length(timepoints$element)){ 
    if(length(sample[[5]][[timepoints$element[i]]][[1]][sample[[5]][[timepoints$element[i]]][[1]]>mzrange[1] & sample[[5]][[timepoints$element[i]]][[1]]<mzrange[2]])!=0){       
      matrixscan<-data.frame(mz=sample[[5]][[timepoints$element[i]]][[1]][sample[[5]][[timepoints$element[i]]][[1]]>mzrange[1] & sample[[5]][[timepoints$element[i]]][[1]]<mzrange[2]],
                             int=sample[[5]][[timepoints$element[i]]][[2]][sample[[5]][[timepoints$element[i]]][[1]]>mzrange[1] & sample[[5]][[timepoints$element[i]]][[1]]<mzrange[2]],
                             timeofscan=timepoints$timeofscan[i])
      if(length(matrixscan$mz)==1) eic[i,]<-matrixscan
      if(length(matrixscan$mz)==0) eic[i,]<-data.frame(mz=paste(mean(mzrange)),int=0,timeofscan=timepoints$timeofscan[i])
      if(length(matrixscan$mz)>1){  
        matrixscan<-matrixscan[matrixscan$int>0,]
        matrixscan<-matrixscan[which.min(abs(matrixscan$mz-mean(mzrange))),]
        eic[i,]<-matrixscan
      }
    } else {
      eic[i,]<-eic[i,]<-data.frame(mz=0,int=0,timeofscan=timepoints$timeofscan[i])
    }
  }
  
  if(ploteic==TRUE){
    if(timerange_old[1]==timerange[1] & timerange_old[2]==timerange[2]){
    title<-paste(extra_title[1],"mz =", format(mean(mzrange),4),"\u00B1",paste(accuracy) ,extra_title[2])
    plot(eic$timeofscan,eic$int,
         xlab="Retention Time (min)",
         ylab="Intensity",
         main=title,...)
    axis(1, at=seq(0,max(eic$timeofscan), by=1))
    } else {
      title<-paste(extra_title[1],"mz =", format(mean(mzrange),4),"\u00B1",paste(accuracy),"&", "tR =", round(timerange[1],2),",",round(timerange[2],2),extra_title[2])
      plot(eic$timeofscan,eic$int,
           xlab="Retention Time (min)",
           ylab="Intensity",
           main=title,...)
    }
  }
  
  return(eic)  
}