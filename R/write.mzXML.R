#'Returns mzXML list object in mzXML files
#'
#'Converts the mzXML list object to a mzXML file and stores it on the hard drive
#'@param mzXML mzXML list created from read.mzXML
#'@param filename The directory in the hard drive that the mzXML files will be stored
#'@param precision option between 32 bit or 64 bit precision
#'@return Returns an mzXML file back to hard drive
#'
#'@author 
#'Originally this function was part of caMassClass R-package by Jarek Tuszynski 
#'@author
#'Codes maintained by Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
#'@export 

write.mzXML<-function(mzXML, filename, precision=c('32', '64'))
{
  Paste  = function(...) paste(..., sep="", collapse="")
  
  fprintf = function(fp, level, ..., append=TRUE)
  { # helper function
    x = paste(..., sep="")
    if (length(x)==0 || is.null(x)) return(NULL)
    spaces = if (level>0) Paste(rep("  ", level)) else ""
    x = gsub("'", "\"", x)
    cat(spaces, x, file=fp, sep="")
    NULL
  }  # done with local functions
  
  library(XML)
  library(digest)
  precision = match.arg(precision)
  if (!is.character(filename)) stop("read.mzXML: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string
  fp  = file(filename, "w")

  if (is.null(mzXML))
    stop("write.mzXML: Variable mzXML has to be an instance of class mzXML");
  if (attr(mzXML, "class")!="mzXML")
    stop("write.mzXML: Variable mzXML has to be an instance of class mzXML");
  
  #------------------------------------
  # Fill-in required sections if empty
  #------------------------------------
  if (is.null(mzXML$header)) {
    Str = "http://sashimi.sourceforge.net/schema_revision/mzXML_2.1"
    mzXML$header = Paste( "<mzXML xmlns='",Str,"'\n",
                          "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n",
                          "xsi:schemaLocation='", Str, " ", Str, "/mzXML_idx_2.1.xsd'>\n")
  }
  if (is.null(mzXML$parentFile)) {
    mzXML$parentFile = Paste( "    <parentFile filename='file://unknown' ",
                              "fileType='RAWData' fileSha1='0000000000000000000000000000000000000000'/>\n")
  }
  if (is.null(mzXML$dataProcessing)) {
    Version = packageDescription("caMassClass")$Version
    Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    mzXML$dataProcessing = Paste("    <dataProcessing>\n",
                                 "      <software type='processing' name='cran.r-project.org/caMassClass' ",
                                 "version='",Version,"' completionTime='",Time,"'/>\n    </dataProcessing>")
  }
  
  #-----------------------------
  # Write beggining of file
  #-----------------------------
  fprintf(fp, 0, "<?xml version='1.0' encoding='ISO-8859-1'?>\n", append=FALSE)
  fprintf(fp, 0, mzXML$header)
  fprintf(fp, 0, "<msRun scanCount='",length(mzXML$scan),"'>\n")
  fprintf(fp, 0, mzXML$parentFile)
  fprintf(fp, 0, mzXML$msInstrument)
  fprintf(fp, 0, mzXML$dataProcessing)
  fprintf(fp, 0, mzXML$separation)
  fprintf(fp, 0, mzXML$spotting)
  
  #---------------------------------
  # Write scan Section
  #---------------------------------
  indexScan = "  <index name='scan'>\n"
  n   = length(mzXML$scan)
  Num = integer(n)
  
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$num
  mzXML$scan = mzXML$scan[ order(Num) ]
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$parentNum
  mzXML$scan = mzXML$scan[ order(Num) ]
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$msLevel
  Num = 1-diff(c(Num,1)) # number of </scan> after each scan

  size = (if (precision=="32") 4 else 8)

  for(i in 1:n) if (!is.null(mzXML$scan[[i]])) {
    indexScan = Paste(indexScan, "    <offset id='",mzXML$scan[[i]]$num,"'>",
                      seek(fp),"</offset>\n")
    mass  = mzXML$scan[[i]]$mass
    peaks = mzXML$scan[[i]]$peaks
    stopifnot(length(mass)==length(peaks))
    if(length(mass)==0){
      peaks<-0
      mass<-0
    }
    
      ioncur<-strsplit(strsplit(mzXML$scan[[i]]$scanAttr, "totIonCurrent=[\"]")[[1]][2], "[\"]")[[1]][1]
      ioncur<-as.numeric(ioncur)
      
      collision_energy_value<-strsplit(strsplit(mzXML$scan[[i]]$scanAttr, "collisionEnergy=[\"]")[[1]][2], "[\"]")[[1]][1]
      if(is.na(collision_energy_value)) collision_energy_value<-25
      
      ScanHeader = Paste("<scan num='"   , mzXML$scan[[i]]$num,
                         "' retentionTime='"   , strsplit(strsplit(mzXML$scan[[i]]$scanAttr, "retentionTime=[\"]")[[1]][2],"[\"]")[[1]][1], 
                         "' polarity='"   , strsplit(strsplit(mzXML$scan[[i]]$scanAttr, "polarity=[\"]")[[1]][2], "[\"]")[[1]][1],
                         "' msLevel='"   , mzXML$scan[[i]]$msLevel, 
                         "' collisionEnergy='"   , collision_energy_value,
                         "' peaksCount='", length(peaks),
                         "' lowMz='"     , sprintf("%.8f",min(mass)),
                         "' highMz='"    , sprintf("%.8f",max(mass)),
                         "' basePeakMz='"    , sprintf("%.8f",mass[which.max(peaks)]),
                         "' basePeakIntensity='"    , sprintf("%.8f",max(peaks)),
                         "' totIonCurrent='"    , sprintf("%.8f",ioncur),
                         "' D'"    , 1)
      
      
      
      
      
      ScanHeader<-substr(ScanHeader,1,nchar(ScanHeader)-4)
      ScanHeader<-Paste(ScanHeader, ">\n")
      
      fprintf(fp, 0, ScanHeader)
      #fprintf(fp, 3, mzXML$scan[[i]]$scanOrigin)
      fprintf(fp, 0, mzXML$scan[[i]]$precursorMz)
      fprintf(fp, 0, mzXML$scan[[i]]$maldi)
      p = as.vector(rbind(mass,peaks))
      fprintf(fp, 0, Paste("<peaks precision='",precision,
                           "' byteOrder='network' pairOrder='m/z-int'>",
                           base64encode(p, endian="big", size=size), "</peaks>\n"))
      fprintf(fp, 0, mzXML$scan[[i]]$nameValue)
      if(Num[i]) for(j in 1:Num[i]) fprintf(fp, 0, "</scan>\n")
    }
#    fprintf(fp, 0, "</scan>\n")
    indexScan = Paste(indexScan, "  </index>\n")
  
  
  #---------------------------------
  # Write end of file
  #---------------------------------
  fprintf(fp, 0, "</msRun>\n")
  n = seek(fp)
  fprintf(fp, 0, indexScan)
  fprintf(fp, 0, "<indexOffset>",n,"</indexOffset>\n")
  cat("  <sha1>", file=fp, sep="")
  n = seek(fp)
  close(fp)
  sha1 = digest(filename, algo="sha1", file=TRUE, length=n)
  cat(sha1, "</sha1>\n</mzXML>\n", file=filename, append=TRUE, sep="")
  invisible(NULL)
  }