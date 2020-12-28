#'Calibrates mzXML list object based on the calibrant peaks which is injected in the beginning (or at the end of) each chromatographic run.
#'
#'Calibrates mzXML list object and returns a calibrated mzXML list on the global environment variables. Then this calibrated mzXML list can be stored as mzXML file using function write.mzXML
#'@param mzXML sample mzXML list created from read.mzXML function
#'@param calibration_region The retention time range in minutes, in which calibrant peaks appear
#'@param calibration_substance It's the calibrant substance, based on which reference mass peaks are calculated.
#'@param int_thres Intensity threshold above which calibrant peaks are taken into account
#'@param mzthres Threshold indicating the absolute difference of experimental and theoretical mass. If more than this value, the calibrant peak is disregarded.
#'@details This function isolates the calibrant peaks, which appear in the beginning of each chromatographic run. Calculates
#'the absolute error in mass Based on the reference masses of the calibrant. After that, an optimized polyonimal model is fitted
#'to the values (error=f(masses)). Based on this expression, all the masses in all the full scans are replaced with the predicted
#'values. 
#'At the moment calibration_substance argument can accept 16 calibrant substances;
#' \itemize{
#'  \item{"Arginine_neg"}
#'  \item{"Arginine_pos"}
#'  \item{"Cs_Perfluoroheptanoate_pos"}
#'  \item{"Fatty_Acids_neg"}
#'  \item{"Fatty_Acids_pos"}
#'  \item{"Li_Formate_neg"}
#'  \item{"Li_Formate_pos"}
#'  \item{"Na_Acetate_neg"}
#'  \item{"Na_Acetate_pos"}
#'  \item{"Na_Formate_neg"}
#'  \item{"Na_Formate_pos"}
#'  \item{"Na_TFA_neg"}
#'  \item{"Na_TFA_pos"}
#'  \item{"PEG_pos"}
#'  \item{"Reserpine_pos"}
#'  \item{"silicone_pos"}
#' }
#'
#'In case another substance is used as calibrant, calibrant peaks can be inserted manually. For example the following
#'format should be followed, so that calibrant peaks can inserted manually:
#'
#'calibration_substance=data.frame(substance=c("Peak1", "Peak2", "Peak3","Peak4","Peak5","Peak6","Peak7"),
#'                                 refmz=c(90.9766,158.9641,226.9515,294.9389,362.9263,430.9138,498.9012))
#'@return calibrated mzXML list object in the global environment
#'
#'@examples
#'mzXMLexample<-read.mzXML(paste(list.files(file.path(find.package("peakTrams"), "data"),pattern=".mzXML",full.names = T)))
#'mzXMLexample_calibrated<-calibration(mzXML=mzXMLexample,calibration_region=c(0.08,0.26), calibrant_substance="Na_Formate_pos",int_thres=1000, mzthres=0.02)
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 
calibration<-function(mzXML=read.mzXML("C:/mymzXMLfolder/sample1.mzXML"),calibration_region=c(0.08,0.26), calibrant_substance="Na_Formate_pos",int_thres=1000, mzthres=0.05){

calibration_region_sec<-calibration_region*60


data("calibrants")
if(any(calibrant_substance==names(calibrants))){
  reference_peaks<-as.data.frame(calibrants[which(calibrant_substance==names(calibrants))])
  names(reference_peaks)<-c("substance","refmz")
} else { 
   reference_peaks<-calibrant_substance 
   }


info<-getinfo(mzXML)
restricted_info<-info[info$timeofscan<calibration_region_sec[2] & info$timeofscan>calibration_region_sec[1],]
restricted_info<-restricted_info[restricted_info$mslevel==1,]

i<-1
j<-1
scanlist<-list()
whichrow<-c(rep(NA,length(reference_peaks$refmz)))

for(i in 1:length(restricted_info$scan)){ #repeat for 5 scans
  spectrum<-data.frame(mz=mzXML$scan[[restricted_info$scan[i]]][[1]],int=mzXML$scan[[restricted_info$scan[i]]][[2]])
  spectrum<-spectrum[spectrum$int>int_thres,]
  
  for(j in 1:length(reference_peaks$refmz)){
  if(min(abs(spectrum[,1]-reference_peaks$refmz[j]))<0.1) whichrow[j]<-which.min(abs(spectrum[,1]-reference_peaks$refmz[j]))
  }
  scanlist[[i]]<-spectrum[c(whichrow),]
}

k<-1
for(k in 1:length(scanlist)) scanlist[[k]][is.na(scanlist[[k]])]<-0

k<-1
sumscanlist<-matrix(nrow=c(length(scanlist)+1),ncol=length(reference_peaks$refmz))
for(k in 1:length(scanlist)) sumscanlist[k,]<-scanlist[[k]][,2]

sumscanlist[c(length(restricted_info[,1])+1),]<-colSums(sumscanlist, na.rm = TRUE, dims = 1)

k<-1
for(k in 1:length(sumscanlist[,1])) sumscanlist[k,]<-sumscanlist[k,]/sumscanlist[length(sumscanlist[,1]),]




k<-1
mzscanlist<-matrix(nrow=c(length(scanlist)+1),ncol=length(reference_peaks$refmz))
for(k in 1:length(scanlist)) mzscanlist[k,]<-scanlist[[k]][,1]

k<-1
avgmass<-c(rep(0,length(mzscanlist[1,])))
for(k in 1 :length(mzscanlist[1,])) avgmass[k]<-sum(mzscanlist[,k]*sumscanlist[,k],na.rm=T)

avgmass[avgmass==0]<-NA
reference_peaks2<-reference_peaks$refmz[!is.na(avgmass)]
avgmass<-avgmass[!is.na(avgmass)]

x=avgmass

i<-1
u<-1
eliminate<-c()
for(i in 1:length(x)){
if(min(abs(x[i]-reference_peaks2))>mzthres){
  eliminate[u]<-i
  u<-u+1
} 
}
if(length(eliminate)!=0){
x<-x[-eliminate]
reference_peaks2<-reference_peaks2[-eliminate]
}

#if(length(x2)!=length(reference_peaks2)){
#reference_peaks2<-reference_peaks2[is.element(round(reference_peaks2),round(x2,0))]
#}

k<-1
y<-c(rep(NA,length(x)))
for(k in 1:length(x)) y[k]<-c(((reference_peaks2[k]-x[k])/min(reference_peaks2[k],x[k])))
y

plot(x,y,main=c("Calibration Curve"),
     ylab="Absolute Error", xlab="m/z", type="p", pch=1, col="red", cex=1, lwd=6)

polyfit <- function(i) x <- AIC(lm(y~poly(x,i)))
opt_order<-as.integer(optimize(polyfit,interval = c(1,length(x)-1))$minimum)
fit1<-lm(y ~ poly(x, opt_order, raw=TRUE))
lines(x, predict(fit1, data.frame(x=x)), col='blue', lwd=2)
###


calibration_info<-info[!(info$timeofscan<calibration_region_sec[2] & info$timeofscan>calibration_region_sec[1]),]

i<-1
for(i in 1:length(calibration_info[,1])){ 
  tobepredicted_mz<-mzXML$scan[[calibration_info$scan[i]]][[1]]
  corfactor<-as.vector(predict(fit1, data.frame(x=tobepredicted_mz)))
  mzXML$scan[[calibration_info$scan[i]]][[1]]<-c(tobepredicted_mz+tobepredicted_mz*corfactor)
  #print(i)
}
return(mzXML)
}