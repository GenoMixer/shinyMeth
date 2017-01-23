# Preprocessing of methylation data

# required packages for differential methylation analysis
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(FlowSorted.Blood.450k)
library(FlowSorted.CordBlood.450k)
library(FlowSorted.DLPFC.450k)
library(RColorBrewer)
library(wateRmelon)
library(MethylAid)
library(shinyMethyl)
library(car)
library(sva)
library(limma)

# gets the 450k annotation data
anno450k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# sets targets directory and reads in targets/ sample sheet for the experiment
setwd("/media/sugi/NGS_Core_Bonn_1/Extremegroups/data")
projectDirectory <- "/media/sugi/NGS_Core_Bonn_1/Extremegroups/data"
targets <- read.metharray.sheet(projectDirectory, pattern="SampleSheet.csv")

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)

# give the data sample names
sampleNames(rgSet) <- targets$Sample_ID

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# normalize the data using quantile normalization
mSetSq <- preprocessQuantile(rgSet)

# normalize the data using Funnorm
mSetFun <- preprocessFunnorm(rgSet)


# effect of different normalisation strategies
par(mfrow = c(2, 2))
# MDS for preprocessRaw, processQuantile, preprocessFunnorm
mdsPlot(getM(mSetRaw),numPositions=1000,sampGroups=targets$subgroup,main="preprocessRaw",pch=20)
mdsPlot(getM(mSetSq),numPositions=1000,sampGroups=targets$subgroup,main="preprocessQuantile",pch=20 )
mdsPlot(getM(mSetFun),numPositions=1000,sampGroups=targets$subgroup,main="preprocessFunnorm",pch=20)


# data before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet,sampGroups=targets$subgroup,main="Raw", legend=FALSE)
legend("topleft", legend=levels(factor(targets$subgroup)),text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq),sampGroups=targets$subgroup,main="Normalized", legend=FALSE)
legend("topleft", legend=levels(factor(targets$subgroup)),text.col=brewer.pal(8,"Dark2"))

