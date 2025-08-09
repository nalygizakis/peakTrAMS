# peakTrAMS R Package

`peakTrAMS` is an R package designed for peak detection, subtraction, and comparison of mass spectrometry data from HILIC and RP chromatographic methods. This README provides a step-by-step guide for running a typical workflow using the `peakTrAMS` package.

---

## Installation

To use the `peakTrAMS` package, install it as follows:

```r
# Install from CRAN or GitHub (example, replace with actual source)
# install.packages("peakTrAMS")
library(peakTrAMS)
```

---

## Workflow Overview

### 1. Specify Paths for Raw Files
Define the paths to the raw `.mzXML` files for HILIC and RP samples and blanks:
```r
path_HILIC <- list.files(paste(find.package(package = "peakTrAMS"), "data", sep = "/"), 
                         pattern = c("Sample_HILIC", ".mzXML"), full.names = TRUE)

path_RP <- list.files(paste(find.package(package = "peakTrAMS"), "data", sep = "/"), 
                      pattern = c("Sample_RP", ".mzXML"), full.names = TRUE)

path_Blank_HILIC <- list.files(paste(find.package(package = "peakTrAMS"), "data", sep = "/"), 
                               pattern = c("Blank_HILIC", ".mzXML"), full.names = TRUE)

path_Blank_RP <- list.files(paste(find.package(package = "peakTrAMS"), "data", sep = "/"), 
                            pattern = c("Blank_RP", ".mzXML"), full.names = TRUE)
```

---

### 2. Load Library
```r
library("peakTrAMS")
```

---

### 3. Read, Calibrate mzXML, and Restrict Processing Area
Create an output directory and calibrate the raw files:
```r
dir.create("output", showWarnings = TRUE)
setwd(file.path(getwd(), "output"))

sample_RP <- calibration(mzXML = read.mzXML(filename = path_RP),
                         calibration_region = c(0.08, 0.26), 
                         calibrant_substance = "Na_Formate_pos", 
                         int_thres = 1000)

blank_RP <- calibration(read.mzXML(filename = path_Blank_RP),
                        calibration_region = c(0.08, 0.26), 
                        calibrant_substance = "Na_Formate_pos", 
                        int_thres = 1000)

sample_HILIC <- calibration(read.mzXML(path_HILIC), 
                            calibrant_substance = "Na_Formate_pos",
                            calibration_region = c(0.08, 0.26), 
                            int_thres = 1000)

blank_HILIC <- calibration(mzXML = read.mzXML(path_Blank_HILIC),
                           calibration_region = c(0.08, 0.26), 
                           calibrant_substance = "Na_Formate_pos", 
                           int_thres = 1000)
```

Restrict the scans based on retention time:
```r
sample_HILIC <- removescans(sample_HILIC, scansORtime = c("beginning", 0.26), time = TRUE)
blank_HILIC <- removescans(blank_HILIC, scansORtime = c("beginning", 0.26), time = TRUE)

sample_RP <- removescans(sample_RP, scansORtime = c("beginning", 0.26), time = TRUE)
blank_RP <- removescans(blank_RP, scansORtime = c("beginning", 0.26), time = TRUE)
```

---

### 4. Scan-by-Scan Subtraction
Subtract blank samples from the corresponding raw samples:
```r
sub_HILIC <- subtract(sample = sample_HILIC, blank = blank_HILIC, mzdiff = 0.01, filter = 100)
sub_RP <- subtract(sample = sample_RP, blank = blank_RP, mzdiff = 0.01, filter = 100)
```

---

### 5. Write Subtracted Chromatograms
Save the subtracted chromatograms to `.mzXML` files:
```r
write.mzXML(sub_RP, "Subtracted_RP.mzXML", precision = "64")
write.mzXML(sub_HILIC, "Subtracted_HILIC.mzXML", precision = "64")
```

---

### 6. Peak Picking and Comparing Peak Lists
Perform peak picking and compare peaks between HILIC and RP:
```r
output <- compareRPHILIC(HILIC = paste(getwd(), "Subtracted_HILIC.mzXML", sep = "/"),
                         RP = paste(getwd(), "Subtracted_RP.mzXML", sep = "/"), 
                         mzthr = 0.01, doubleckeck = TRUE, 
                         ppm = 30, peakwidth = c(14.34, 50), sn = 10)
```

---

### 7. Prioritization of Common and Uncommon Peaks
Prioritize common and unique peaks:
```r
prioritized_output <- prioritization(output)
```

---

### 8. Plot Extracted Ion Chromatograms (EICs)
Generate plots for prioritized peaks:
```r
writeoutput(prioritized_output, sampleRP = sub_RP, sampleHILIC = sub_HILIC, type = "o")
```

This will create PDF files containing the EICs of unique and common peaks.

---

## Outputs
- `Subtracted_RP.mzXML` and `Subtracted_HILIC.mzXML`: Subtracted chromatograms.
- PDFs with unique peaks in HILIC and RP as well as prioritized common peaks.

---

## Notes
- Ensure the `.mzXML` files are formatted correctly.
- Adjust parameters (`mzdiff`, `ppm`, `peakwidth`, etc.) as needed for your dataset.

---
