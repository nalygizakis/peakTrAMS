#'Subtracts procedural blank from sample chromatogram
#'
#'This function takes into account only MS1 levels (other levels are removed). 
#'It goes to all scans of sample, searches for the nearest in retention time full scan in blank and 
#'then reads the spectral peaks. Aftewards, checks one by one the spectral peaks in sample and 
#'searches for possible match in blank withing a given mass accuracy threshold (see argument mzdiff). 
#'In case a spectral peak is present in blank then the intensity of the sample is subtracted by the 
#'intenisity observed in blank. In case the result is a negative intensity value, it is converted to zero.
#'@param sample sample mzXML list (provided from read.mzXML function)
#'@param blank blank mzXML list (provided from read.mzXML function)
#'@param mzdiff Absolute mass accuracy threshold
#'@param filter cut off filter (normally is set at a small intensity value)
#'
#'@return an mzXML list object containg subtracted chromatogram
#'
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 
subtract<-function(sample=read.mzXML(paste(list.files(file.path(find.package("peakTrams"), "data"),pattern=c("Sample_RP",".mzXML"),full.names = T))),
                   blank=read.mzXML(paste(list.files(file.path(find.package("peakTrams"), "data"),pattern=c("Blank_RP",".mzXML"),full.names = T))),
                   mzdiff=0.01, filter=100){

    i<-1
    k<-1
    numsample<-data.frame(element_list=0,scan=0,ret_time=0)
    for(i in 1:length(sample$scan)){
      if(sample[[5]][[i]]$msLevel==1){
      numsample[k,]<-c(i,sample[[5]][[i]]$num,as.numeric(strsplit(strsplit(sample[[5]][[i]]$scanAttr,split="S")[[1]][1],split="PT")[[1]][2]))
      k<-k+1
    }   
  }
     cat("Sample has", paste(length(sample$scan)), "scans from which", paste(k-1),"scans are in MS1 level and", c(length(sample$scan)-k+1),"in MS2 level. MS2 scans are disregarded during the subtraction procedure","\n")
    
    i<-1
    u<-1
    numblank<-data.frame(element_list=0,scan=0,ret_time=0)
    for(i in 1:length(blank$scan)){
      if(blank[[5]][[i]]$msLevel==1){
        numblank[u,]<-c(i,blank[[5]][[i]]$num,as.numeric(strsplit(strsplit(blank[[5]][[i]]$scanAttr,split="S")[[1]][1],split="PT")[[1]][2]))
        u<-u+1
      }   
    }
    
    cat("Blank has", paste(length(blank$scan)), "scans from which", paste(u-1),"scans are in MS1 level and", c(length(blank$scan)-u+1),"in MS2 level. MS2 scans are disregarded during the subtraction procedure","\n")
    
    progress<-txtProgressBar(min=1, max=length(numsample[,1]), style=3)
    
    scan<-1
    for(scan in 1:length(numsample[,1])){
    scan2<-which.min(abs(numsample$ret_time[scan]-numblank$ret_time))
        
    matrix_sample<-cbind(sample[[5]][[numsample$element_list[scan]]][[1]][],sample[[5]][[numsample$element_list[scan]]][[2]][])
    matrix_blank<-cbind(blank[[5]][[numblank$element_list[scan2]]][[1]][],blank[[5]][[numblank$element_list[scan2]]][[2]][]) 
    matrix_subtracted<-matrix(0,nrow=length(matrix_sample[,1]), ncol=2)

    if(length(matrix_sample[,1])!=0 & length(matrix_blank[,1])!=0){
    i<-1
    for(i in 1:length(matrix_sample[,1])){
        if(min(abs(matrix_sample[i,1]-matrix_blank[,1]))<mzdiff){
          temp<-which.min(abs(matrix_sample[i,1]-matrix_blank[,1]))
          matrix_subtracted[i,]<-c(matrix_sample[i,1], matrix_sample[i,2]-matrix_blank[temp,2])
        } else {
          matrix_subtracted[i,]<-c(matrix_sample[i,1],matrix_sample[i,2])
        }
    }

    
    matrix_subtracted<-matrix_subtracted[matrix_subtracted[,2]>filter,]
    ##Check with a plot that the function works
    ##Use scan=309
    #par(mfrow=c(3,1))
    #plot(matrix_sample,type="h",main=numsample$scan[scan],xaxt="n")
    #axis(1, at=seq(from=0,to=1000, by=50))
    #plot(matrix_blank,type="h", main=numblank$scan[scan2],xaxt="n")
    #axis(1, at=seq(from=0,to=1000, by=50))
    #plot(matrix_subtracted,type="h",xaxt="n")
    #axis(1, at=seq(from=0,to=1000, by=50))
    
    if(!is.vector(matrix_subtracted)){
    sample[[5]][[numsample$element_list[scan]]][[1]]<-matrix_subtracted[,1]
    sample[[5]][[numsample$element_list[scan]]][[2]]<-matrix_subtracted[,2]
    } else {
      sample[[5]][[numsample$element_list[scan]]][[1]]<-matrix_subtracted[1]
      sample[[5]][[numsample$element_list[scan]]][[2]]<-matrix_subtracted[2]
    }
    

    }
    setTxtProgressBar(progress, scan)
    }
      return(sample)
}