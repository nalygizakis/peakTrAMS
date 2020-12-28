#'Performs Extracted Ion Chromatograms (EICs) for unique peaks in RP, unique peaks in HILIC, and common peaks of category 5.
#'
#'Plots the prioritized output and saves in pdf files the EICs. User should specify the path of the mzXML files. In case numofcommon,
#'numofuniqueHILIC and numofuniqueRP are not passed then all the EICs will be plotted, which may lead to high waiting time.
#'@param prioritized_output output created by prioritization function
#'@param sampleRP mzXML list object of chromatograms in RP
#'@param sampleHILIC mzXML list object of chromatograms in HILIC
#'@param accuracy absolute mass accuracy
#'@param numofcommon number of common features (from category 5) that user wants to be plotted
#'@param numofuniqueHILIC number of HILIC unique features that user wants to be plotted
#'@param numofuniqueRP number of RP unique features that user wants to be plotted 
#'
#'@return Created 3 csv files with common peaks, unique peaks in HILIC and 
#'unique peaks in RP and 3 pdf files with the EICs.
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 
writeoutput<-function(prioritized_output=prioritized_output,
                      sampleRP=paste(getwd(),"Subtracted_RP.mzXML",sep="/"),
                      sampleHILIC=paste(getwd(),"Subtracted_HILIC.mzXML",sep="/"),
                      accuracy=0.005,
                      numofcommon=length(prioritized_output[[1]][prioritized_output[[1]]$category==5,1]),
                      numofuniqueHILIC=length(prioritized_output[[2]][,1]),
                      numofuniqueRP=length(prioritized_output[[3]][,1]),...){
    
  
 
   # sampleRP<-read.mzXML(paste(getwd(),"Subtracted_HILIC.mzXML",sep="/"))
   # sampleHILIC<-read.mzXML(paste(getwd(),"Subtracted_RP.mzXML",sep="/"))
    
    write.csv(prioritized_output[[1]],"common.csv")
    write.csv(prioritized_output[[2]],"unique_HILIC.csv")
    write.csv(prioritized_output[[3]],"unique_RP.csv")
    
    
    ########################################################################
    pdf("uniquepeaksHILIC.pdf")
    i<-1
    cat("Creating pdf file with", numofuniqueHILIC, "EICs of unique peaks in HILIC","\n")
    progress<-txtProgressBar(min=1, max=numofuniqueHILIC, style=3)
    for(i in 1:numofuniqueHILIC){
      par(mfrow=c(2,1))
      plotXIC(sample=sampleRP,mzrange=prioritized_output[[2]][i,1],accuracy=accuracy,
                         ...,extra_title = "RP, EIC of")
      eicsub_hilic<-plotXIC(sampleHILIC,mzrange=prioritized_output[[2]][i,1],accuracy=accuracy,
                            ...,extra_title = "HILIC, EIC of")
      rect(prioritized_output[[2]][i,4]/60-0.5, 0, prioritized_output[[2]][i,4]/60+0.5, max(eicsub_hilic), density = NULL, angle = 45)
      setTxtProgressBar(progress, i)
    }
    cat("\n","Waiting...Writing pdf file","\n")
    dev.off()
    
    
    ########################################################################
    pdf("uniquepeaksRP.pdf")
    i<-1
    cat("Creating pdf file with", numofuniqueRP, "EICs of unique peaks in RP","\n")
    progress<-txtProgressBar(min=1, max=numofuniqueRP, style=3)
    for(i in 1:numofuniqueRP){
      par(mfrow=c(2,1))
      eicsub_RP<-plotXIC(sampleRP,mzrange=prioritized_output[[3]][i,1],accuracy=accuracy,
              ...,extra_title = "RP, EIC of")
      rect(prioritized_output[[3]][i,4]/60-0.5, 0, prioritized_output[[3]][i,4]/60+0.5, max(eicsub_RP), density = NULL, angle = 45)
      plotXIC(sampleHILIC,mzrange=prioritized_output[[3]][i,1],accuracy=accuracy,
              ..., extra_title = "HILIC, EIC of")
      setTxtProgressBar(progress, i)
    }
    cat("\n","Waiting...Writing pdf file","\n")
    dev.off()
    ######################################################################
    
    
    pdf("common_peaks_Category5.pdf")
    i<-1
    
    cat("\n","Creating pdf file with", numofcommon, "EICs of common peaks","\n")
    progress2<-txtProgressBar(min=1, max=numofcommon, style=3)
    for(i in 1:c(numofcommon)){
      par(mfrow=c(2,1))
      eicsub_rp<-plotXIC(sampleRP,mzrange=prioritized_output[[1]][i,1],accuracy=accuracy,
                         ...,extra_title = "RP, EIC of")
      rect(prioritized_output[[1]][i,3]/60-0.5, 0, prioritized_output[[1]][i,3]/60+0.5, max(eicsub_rp), density = NULL, angle = 45)
      eicsub_hilic<-plotXIC(sampleHILIC,mzrange=prioritized_output[[1]][i,1],accuracy=accuracy,
              ...,extra_title = "HILIC, EIC of")
      rect(prioritized_output[[1]][i,4]/60-0.5, 0, prioritized_output[[1]][i,4]/60+0.5, max(eicsub_hilic), density = NULL, angle = 45)
      setTxtProgressBar(progress2, i)
    }
    cat("\n","Waiting...Writing pdf file","\n")
    dev.off()
    ######################################################################
    return(NULL)
  }