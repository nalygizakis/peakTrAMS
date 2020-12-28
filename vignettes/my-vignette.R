## ---- fig.show='hold'----------------------------------------------------
  #1.Specify the paths of the raw files
  path_HILIC<-list.files(paste(find.package(package="peakTrAMS"),"data",sep="/"),pattern = c("Sample_HILIC",".mzXML"), full.names = TRUE)
  path_RP<-list.files(paste(find.package(package="peakTrAMS"),"data",sep="/"),pattern = c("Sample_RP",".mzXML"), full.names = TRUE)
  path_Blank_HILIC<-list.files(paste(find.package(package="peakTrAMS"),"data",sep="/"),pattern = c("Blank_HILIC",".mzXML"), full.names = TRUE)
  path_Blank_RP<-list.files(paste(find.package(package="peakTrAMS"),"data",sep="/"),pattern = c("Blank_RP",".mzXML"), full.names = TRUE)

  rpath_HILIC<-"E:/4.Peaktrams/dataset_2015B/In_07.03.2015_AutoMS2_HILIC_GD7_01_12978.mzXML"
  path_RP<-"E:/4.Peaktrams/dataset_2015B/In_07.03.2015_AutoMS2_RP_GD7_01_12638.mzXML"
  path_Blank_HILIC<-"E:/4.Peaktrams/dataset_2015B/Procedural_Blank_HILIC_GB1_01_12893.mzXML"
  path_Blank_RP<-"E:/4.Peaktrams/dataset_2015B/Procedural_Blank_RP_GB1_01_12559.mzXML"
  
    #1.Load libraries needed
  library("peakTrAMS")
  
  #3.Read, calibrate mzXML and restrict processing area
  dir.create(c("output"),showWarnings=T)
  setwd(as.character(paste(getwd(),c("output"), sep="/")))
  
  sample_RP<- calibration(mzXML=read.mzXML(filename=path_RP),calibration_region=c(0.08,0.26), calibrant_substance="Na_Formate_pos",int_thres=1000)
  blank_RP<- calibration(read.mzXML(filename=path_Blank_RP),calibration_region=c(0.08,0.26), calibrant_substance="Na_Formate_pos",int_thres=1000)
  sample_HILIC<- calibration(read.mzXML(path_HILIC),calibration_region=c(0.08,0.26), calibrant_substance="Na_Formate_pos",int_thres=1000)
  blank_HILIC<- calibration(mzXML=read.mzXML(path_Blank_HILIC),calibration_region=c(0.08,0.26), calibrant_substance="Na_Formate_pos",int_thres=1000)
  
  
  blank_HILIC<-removescans(mzXML=blank_HILIC,scansORtime=c("beginning",0.26),time=TRUE)
  blank_HILIC<-removescans(mzXML=blank_HILIC,scansORtime=c(17,"end"),time=TRUE)
  sample_HILIC<-removescans(mzXML=sample_HILIC,scansORtime=c("beginning",0.26),time=TRUE)
  sample_HILIC<-removescans(mzXML=sample_HILIC,scansORtime=c(17,"end"),time=TRUE)
  sample_RP<-removescans(mzXML=sample_RP,scansORtime=c("beginning",0.26),time=TRUE)
  blank_RP<-removescans(mzXML=blank_RP,scansORtime=c("beginning",0.26),time=TRUE)
  
  
  #4.Scan-by-scan Subsraction
  sub_HILIC<-subtract(sample=sample_HILIC,blank=blank_HILIC,mzdiff=0.01,filter=100) 
  sub_RP<-subtract(sample=sample_RP,blank=blank_RP,mzdiff=0.01,filter=100)
  
  
  #5.Write the subtracted chromatograms
  write.mzXML(sub_RP, paste("Subtracted_RP",".mzXML", sep=""), precision=c('64'))
  write.mzXML(sub_HILIC, paste("Subtracted_HILIC",".mzXML", sep=""), precision=c('64'))
  
  #6.Peak picking and comparing of peaklists
  output<-compareRPHILIC(HILIC="C:/Users/Acer/OneDrive/R_working_folder/peakTrams_new8_2016/output/output/Subtracted_HILIC.mzXML",
                 RP="C:/Users/Acer/OneDrive/R_working_folder/peakTrams_new8_2016/output/output/Subtracted_RP.mzXML", 
                 mzthr=0.015, doubleckeck=TRUE,
                 ppm=30, peakwidth=c(30,61.5), sn=10)
  
 
  #7.Prioritization of common and uncommon peaks
  prioritized_output<-prioritization(output)
  
  library("VennDiagram", lib.loc="~/R/win-library/3.2")
  venn.plot <- venn.diagram(
    list(RP = 1:c(length(prioritized_output[[3]][,1])+length(prioritized_output[[1]][,1])*2), 
         HILIC = length(prioritized_output[[3]][,1]):c(length(prioritized_output[[3]][,1])+length(prioritized_output[[1]][,1])*2+length(prioritized_output[[2]][,1]))),
    "mzthr0.01_ppm17.6_width14-50_sn10_HILIC_10Xmzthr+noise300.tiff"
    )

  #8.Plot EICs
  writeoutput(prioritized_output,
              sampleRP=read.mzXML(RP),
              sampleHILIC=read.mzXML(HILIC),
              type="o")
  
  
library(xcms)  
aa<-xcmsSet(HILIC,method = 'centWave',
    ROI.list=list(data.frame(scmax=10000000,scmin=1,
                             mzmin=c(838.5369-0.04),
                             mzmax=c(838.5369+0.04))))
  aa@peaks
  
  