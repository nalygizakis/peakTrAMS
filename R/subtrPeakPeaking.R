#'sample subtraction based on accurate mass, retention time and maximum intensity
#'
#'Due to chromatographic drift scan by scan subtraction may fail to remove completely blank peaks....bla bla bla
#'@param should write this ...
#'@return should write this....
#'@export 
subtrPeakPeaking<-function (sample, blank, mzdiff = 0.01, rtrangediff = 20, times=3){ 
  i<-1
  j<-1
  k<-1
  remove<-list()
  sample2<-as.data.frame(sample@peaks) 
  sample<-as.data.frame(sample@peaks)
  blank<-as.data.frame(blank@peaks)
  sink("features_removed.txt")
  for(i in 1:length(sample$mz)){
    if((min(abs((sample$mz[i] - blank$mz))))<mzdiff){
      
      blank_temp<-blank$mz
      howmany<-length(abs(sample$mz[i]-blank$mz)[c(sort(abs(sample$mz[i]-blank$mz),decreasing=F)<mzdiff)])
      temp<-c()
      temp[j]<-which.min(abs(sample$mz[i]-blank$mz))
      if(howmany>1){ 
        for(j in 1:(howmany-1)){ 
          blank_temp[c(temp[j])]<-10^9
          temp[j+1]<-which.min(abs(sample$mz[i]-blank_temp))
        }
      }
      j<-1



      
      for(k in 1:length(temp)){
        range <- abs(sample$rt[i]-blank$rt[temp[k]])
        if (range < rtrangediff){       
          if(blank$maxo[temp[k]]*times>sample$maxo[i]){ 
          remove[[i]]<-i
          cat("==============================","\n")
          cat("Feature from sample","\n")
          print(sample[c(i),])
          cat("was removed, because in blank there is","\n")
          print(blank[c(temp[k]),])
          } 
        }
      } 
      k<-1

    }
}

remove<-c(do.call(rbind.data.frame, remove))

if(length(remove)==0){
  print("No peaks removed by subtraction of peaklists.")
} else { 
remove<-unique(remove[[1]])
sample2<-sample2[-c(remove),]


    print(sprintf("Initially it was: %s features and remained %s features. Lines removed were %s", 
                  length(sample$mz), length(sample2$mz), (length(sample$mz)-length(sample2$mz))))
    print(sprintf("In detail Features removed were %s", 
                   paste(remove, collapse = ",")))
}
closeAllConnections() 
output<-data.frame(mz=sample2[,c(-4,-2,-3,-5:-11)],
                        rt=round(sample2[,c(-1,-2,-3,-5:-11)]/60,2),
                        maxpeak_area_HILIC=sample2[,c(-1,-4,-2,-3,-5,-6,-7,-8,-10,-11)])



    return(output)
  } 