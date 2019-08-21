# Quality control, probe filtering and sex check

setwd("/media/sugi/NGS_Core_Bonn_1/Extremegroups/data")

# starts interactive visualization using shinyMethyl 
summary <- shinySummarize(rgSet)
summary.norm <- shinySummarize(mSetFun)
runShinyMethyl(summary, summary.norm)

# gets qc plot, removes poor quality samples and generates qc report
qc <- getQC(mSetRaw)
plot1 <- plotQC(qc)
passed_qc = rownames(qc)[(((qc$mMed+qc$uMed)/2)>10.5)]
rgSet = rgSet[,sampleNames(rgSet) %in% passed_qc]
targets = subset(targets, Basename %in% passed_qc)
qcReport(rgSet,sampNames=targets$Sample_ID,sampGroups=targets$Sample_Group,pdf="20161212_qcReport_extreme_groups.pdf")

# sex prediction
predictedSex <- getSex(mSetSq, cutoff = -2)
plotSex(getSex(mSetSq, cutoff = -2))
predSex = data.frame(Basename=row.names(predictedSex),predsex=predictedSex)
nominalSex = data.frame(Sample_ID=targets$Sample_ID,sex=targets$gender)
checkSex = merge(predSex,nominalSex,by="Sample_ID")
table(checkSex[c("predSex","sex")],exclude=NULL)
write.table(checkSex, "20161212_sexcheck_extreme_groups.txt", quote=F, sep="\t")

# mean detection p-value across all samples
detP <- detectionP(rgSet)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$subgroup)], las=2,cex.names=0.8, ylim = c(0,0.002), ylab="Mean detection p-values")
legend("topright", legend=levels(factor(targets$subgroup)), fill=pal,bg="white",cex=0.75)

# removes poor quality samples
keep <- colMeans(detP) < 0.05
mSetFunFlt <- mSetFun[,keep]
mSetSqFlt <- mSetSq[,keep]
targets <- targets[keep,]
detP <- detP[keep,]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetFun) 
table(keep)
mSetFunFlt <- mSetFun[keep,]
mSetSqFlt <- mSetSq[keep,]

# exclude probes with a beadcount < 3
RGSet_ex <- read.metharray.exp(targets=targets,extended=TRUE) # extended set
beadcount <- beadcount(RGSet_ex)
beadcount_NAs <- apply(beadcount, 1, function(z) sum(is.na(z)))
failed_beads_probes <- beadcount[(beadcount_NAs > 2),]
failed_beads_probes <- row.names(failed_beads_probes)
length(failed_beads_probes)
keep <- !(featureNames(mSetFunFlt) %in% failed_beads_probes)
table(keep)
mSetFunFlt <- mSetFunFlt[keep,]
mSetFunFlt


# remove probes on the sex chromosomes
keep <- !(featureNames(mSetFunFlt) %in% annoEPIC$Name[annoEPIC$chr %in% c("chrX","chrY")])
table(keep)
mSetFunFlt <- mSetFunFlt[keep,]
mSetSqFlt <- mSetSq[keep,]

# remove probes with SNPs at CpG site
mSetFunFlt <- dropLociWithSnps(mSetFunFlt)
mSetSqFlt <- dropLociWithSnps(mSetSq)

# exclude cross reactive probes
xReactiveProbes <- read.csv("mmc2_cross_reactive_probes.csv", sep="/", stringsAsFactors=FALSE)
keep <- !(featureNames(mSetFunFlt) %in% xReactiveProbes$Target_ID)
table(keep)
mSetFunFlt <- mSetFunFlt[keep,]
mSetSqFlt <- mSetSq[keep,]

# exclude non specific probes (ch.)
non_specific_Probes <- read.csv("mmc3_non_cpg_sites.csv", sep="/", stringsAsFactors=FALSE)
keep <- !(featureNames(mSetFunFlt) %in% non_specific_Probes$Target_ID)
table(keep)
mSetFunFlt <- mSetFunFlt[keep,]
mSetSqFlt <- mSetSq[keep,]

# extracts M-values for differential methylation analysis
mqn <- getM(mSetSqFlt)
mfn <- getM(mSetFnFlt)

# saves normalized datasets and target file
save(mSetSqFlt,mSetFnFlt,targets,file="media/sugi/NGS_Core_Bonn_1/Extremegroups/data/methylation_filtered_extremegroups.Rdata",compress=F)






