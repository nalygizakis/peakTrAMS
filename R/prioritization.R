#'Prioritizes unique and common peak lists. 
#'
#'Unique peaks are prioritized based on the observed maximum intensity. 
#'Common peaks are prioritized based on three criteria. 
#'The first criterio is that peaks should be retained in HILIC chromatography. 
#'This happens only if compounds elutes after the dead volume of the column, which can be set 
#'from deadvolume argument.
#'The second criterio is is related to the fact that we are looking for polar compounds. Such compounds
#'are eluted relatively in the beginning of Reversed-phase chromatography, whether are eluted later in HILIC
#'chromatography. So, the second requirement is that retention time in HILIC is higher than in RP.
#'The third criterio states that polar compounds should give better response in HILIC than in RP.
#'
#'According to the 3 criteria, common peaks are divided into 8 categories, which were eliminated to 5
#'because 3 categories have no analytical explanation and normally consist of very few common pair peaks.
#'
#'The categories are the followings:
#' \itemize{
#'  \item{"Category 5"}{ All the criteria are valid. Therefore, these are the peaks of interest}
#'  \item{"Category 4"}{ Criteria 1 and 2 are valid. Still some compounds of interest may be in this group}
#'  \item{"Category 3"}{ Criteria 1 and 3 are valid.}
#'  \item{"Category 2"}{ Criterio 1 is valid.}
#'  \item{"Category 1"}{ Criterio 1 is invalid.}
#' }
#'@param output object returned from compare function.
#'@param deadvolume dead volume of HILIC column
#'@return Returns a list with two data frame objects, which contain the prioritized common and unique peaks respectively.
#'@author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export

prioritization<-function(output,deadvolume=c(1.5)){

  prior_output<-list()
  i<-1
for(i in 1:length(output[[1]])){
  common<-data.frame()
  temp_HILIC<-output[[1]][[i]][output[[1]][[i]]$chromatogr=="HILIC",]
  temp_RP<-output[[1]][[i]][output[[1]][[i]]$chromatogr=="RP",]
 
  combinrettimes<-as.data.frame(expand.grid(temp_RP$rt, temp_HILIC$rt))
  common<-data.frame(mz_RP=0,mz_HILIC=0,rt_RP=0,rt_HILIC=0,
                      into_RP=0,into_HILIC=0,intb_RP=0,intb_HILIC=0,
                      maxo_RP=0,maxo_HILIC=0,
                      sn_RP=0,sn_HILIC=0)
  names(combinrettimes)<-c("rt_RP","rt_HILIC")
  k<-1
  for(k in 1:length(combinrettimes[,1])){
    whichRP<-which(combinrettimes[k,1]==temp_RP$rt)
    whichHILIC<-which(combinrettimes[k,2]==temp_HILIC$rt)
    common[k,]<-c(mz_RP=temp_RP[whichRP,]$mz,mz_HILIC=temp_HILIC[whichHILIC,]$mz,rt_RP=temp_RP[whichRP,]$rt,rt_HILIC=temp_HILIC[whichHILIC,]$rt,
                           into_RP=temp_RP[whichRP,]$into,into_HILIC=temp_HILIC[whichHILIC,]$into,intb_RP=temp_RP[whichRP,]$intb,intb_HILIC=temp_HILIC[whichHILIC,]$intb,
                           maxo_RP=temp_RP[whichRP,]$maxo,maxo_HILIC=temp_HILIC[whichHILIC,]$maxo,
                           sn_RP=temp_RP[whichRP,]$sn,sn_HILIC=temp_HILIC[whichHILIC,]$sn)
  }
  prior_output[[i]]<-common
  #print(i)
}
  
  
  #Prioritizing common peaks
  i<-1
  j<-1
  for(i in 1:length(prior_output)){
    prior_output[[i]]$betResponse<-prior_output[[i]]$maxo_HILIC>prior_output[[i]]$maxo_RP
    prior_output[[i]]$retainedHILIC<-prior_output[[i]]$rt_HILIC>deadvolume*60
    prior_output[[i]]$retHILIChigher<-prior_output[[i]]$rt_HILIC>prior_output[[i]]$rt_RP
    prior_output[[i]]$valid<-as.numeric(prior_output[[i]]$betResponse+prior_output[[i]]$retainedHILIC+prior_output[[i]]$retHILIChigher)
    prior_output[[i]]$category<-NA
    for(j in 1:length(prior_output[[i]][,1])){
      if(prior_output[[i]]$valid[j]==3) prior_output[[i]]$category[j]<-5
      if(prior_output[[i]]$valid[j]==2 & prior_output[[i]]$betResponse[j]==FALSE) prior_output[[i]]$category[j]<-4
      if(prior_output[[i]]$valid[j]==2 & prior_output[[i]]$retHILIChigher[j]==FALSE) prior_output[[i]]$category[j]<-3
      if(prior_output[[i]]$valid[j]==2 & prior_output[[i]]$retainedHILIC[j]==FALSE) prior_output[[i]]$category[j]<-1
      if(prior_output[[i]]$valid[j]==1 & prior_output[[i]]$retHILIChigher[j]==FALSE  & prior_output[[i]]$betResponse[j]==FALSE) prior_output[[i]]$category[j]<-2
      if(prior_output[[i]]$valid[j]==1 & prior_output[[i]]$retainedHILIC[j]==FALSE  & prior_output[[i]]$betResponse[j]==FALSE) prior_output[[i]]$category[j]<-1
      if(prior_output[[i]]$valid[j]==1 & prior_output[[i]]$retainedHILIC[j]==FALSE  & prior_output[[i]]$retHILIChigher[j]==FALSE) prior_output[[i]]$category[j]<-1
      if(prior_output[[i]]$valid[j]==0) prior_output[[i]]$category[j]<-1
            }
    prior_output[[i]]$betResponse<-prior_output[[i]]$maxo_HILIC/prior_output[[i]]$maxo_RP
  }
  
  prior_common<-do.call(rbind.data.frame, prior_output)
  prior_common$mz_RP<-round(prior_common$mz_RP,4)
  prior_common$mz_HILIC<-round(prior_common$mz_HILIC,4)
  prior_common$maxo_RP<-round(prior_common$maxo_RP,0)
  prior_common$maxo_HILIC<-round(prior_common$maxo_HILIC,0)

  prior_common<-prior_common[order(prior_common$valid, decreasing = T),]

  plot.new()
  plot(x=prior_common[,1],
       y=log(c(prior_common[,4]/prior_common[,3])),
       ylab="log(Ret. Time HILIC/Ret. Time RP)",
       xlab="mz",
       main="Common Features",
       pch=prior_common$category,
       col=prior_common$category,
#       xlim=c(80,900),     
       #xlim=c(min(prior_common[,1]),max(prior_common[,1])),
       xaxt = "n")
  axis(1, at=seq(0,max(prior_common[,1]), by=10))
  abline(h = 0, v = 0, col = "gray60")
  
  legend(x="topright", legend=paste("Category",sort(unique(prior_common$category), decreasing = T)), col=sort(unique(prior_common$category), decreasing = T),pch=sort(unique(prior_common$category), decreasing = T))
  
  png("Common_features.png")
  plot.new()
  plot(x=prior_common[,1],
       y=log(c(prior_common[,3]/prior_common[,4])),
       ylab="log(Ret. Time HILIC/Ret. Time RP)",
       xlab="mz",
       main="Common Features",
       pch=prior_common$category,
       col=prior_common$category,
       xlim=c(80,900),     
       #xlim=c(min(prior_common[,1]),max(prior_common[,1])),
       xaxt = "n")
  axis(1, at=seq(0,max(prior_common[,1]), by=10))
  abline(h = 0, v = 0, col = "gray60")
  
  legend(x="topright", legend=paste("Category",sort(unique(prior_common$category), decreasing = T)), col=sort(unique(prior_common$category), decreasing = T),pch=sort(unique(prior_common$category), decreasing = T))
  dev.off()
  
  
  cat("In total",length(prior_common[,1]),"common features distributed in 5 caregories","\n")
  cat("Category 5:",length(prior_common$category[prior_common$category==5]),"\n")
      cat("Category 4:",length(prior_common$category[prior_common$category==4]),"\n")
          cat("Category 3:",length(prior_common$category[prior_common$category==3]),"\n")
              cat("Category 2:",length(prior_common$category[prior_common$category==2]),"\n")
                  cat("Category 1:",length(prior_common$category[prior_common$category==1]),"\n")
                
                  cat("Unique peaks in HILIC",length(output[[2]]),"\n")
                  cat("Unique peaks in RP",length(output[[3]]),"\n")
                  #sink("common_peaks_report.txt")
                  #cat("Unique peaks in HILIC",length(output[[2]]),"\n")
                  #cat("Unique peaks in RP",length(output[[3]]),"\n")
                  #cat("In total",length(prior_common[,1]),"common features distributed in 5 caregories","\n")
                  #cat("Category 5:",length(prior_common$category[prior_common$category==5]),"\n")
                  #cat("Category 4:",length(prior_common$category[prior_common$category==4]),"\n")
                  #cat("Category 3:",length(prior_common$category[prior_common$category==3]),"\n")
                  #cat("Category 2:",length(prior_common$category[prior_common$category==2]),"\n")
                  #cat("Category 1:",length(prior_common$category[prior_common$category==1]),"\n")
                  
                  #closeAllConnections() 
                  
  #Prioritize uncommon peaks
  uniqueHILIC<-do.call(rbind.data.frame,output[[2]])
  uniqueHILIC<-uniqueHILIC[order(uniqueHILIC$maxo ,decreasing = T),]

  uniqueRP<-do.call(rbind.data.frame,output[[3]])
  uniqueRP<-uniqueRP[order(uniqueRP$maxo ,decreasing = T),]

  prioritized_output<-list()
  prioritized_output[[1]]<-prior_common
  prioritized_output[[2]]<-uniqueHILIC
  prioritized_output[[3]]<-uniqueRP

  return(prioritized_output)
}