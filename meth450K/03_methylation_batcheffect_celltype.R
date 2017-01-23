# Batch effect and cell type correction

setwd("/media/sugi/NGS_Core_Bonn_1/Extremegroups/data")

# clean M-Values
varfn = apply(mfn,1,var)
varfq = apply(mqn,1,var)
mfn = mfn[!is.na(varfn),]
mqn = mqn[!is.na(varfq),]

# model design for Combat batch effect correction
mod0 = model.matrix(~1, data=targets)

# batch variables
targets$Date_of_typing <- as.factor(targets$Date_of_typing)
targets$Slide <- as.factor(targets$Slide)
targets$Sentrix_slide <- as.factor(targets$Array)

# iterative batch effect correction with known batch variables
mfn1 = ComBat(mfn,batch=targets$Date_of_typing,mod0)
mqn1 = ComBat(mqn,batch=targets$Date_of_typing,mod0)

mfn2 = ComBat(mfn1,batch=targets$Slide,mod0)
mqn2 = ComBat(mqn1,batch=targets$Slide,mod0)

mfn3 = ComBat(mfn2,batch=targets$Array,mod0)
mqn3 = ComBat(mqn2,batch=targets$Array,mod0)

# number of surrogate variables
num.sv(mfn,mod0,vfilter=10000)
num.sv(mqn,mod0,vfilter=10000)

# executes PCA for most variable sites and extracts scores on first components
pcv = function(x,cols=5,sites=10000) {
  var1 = apply(x,1,var)
  set1 = order(var1,decreasing=TRUE)[1:sites]
  pcs = prcomp(t(x[set1,]))
  return(pcs$x[,1:cols])
}

# control with anova
pcfn1 = pcv(mfn1)
pcqn1 = pcv(mqn1)

pcfn1 = pcv(mfn1)
pcqn1 = pcv(mqn1)

pcfn1 = pcv(mfn1)
pcqn1 = pcv(mqn1)

Anova(lm(pcfn1~factor(targets$Date_of_typing)))
Anova(lm(pcqn1~factor(targets$Date_of_typing)))

Anova(lm(pcfn2~factor(targets$Slide)))
Anova(lm(pcqn2~factor(targets$Slide)))

Anova(lm(pcfn3~factor(targets$Array)))
Anova(lm(pcqn3~factor(targets$Array)))

pairs(pcfn1,col=factor(targets$Date_of_typing),pch=20)
pairs(pcqn1,col=factor(targets$Date_of_typing)),pch=20)

pairs(pcfn2,col=factor(targets$Slide),pch=20)
pairs(pcqn2,col=factor(targets$Slide)),pch=20)

pairs(pcfn3,col=factor(targets$Array),pch=20)
pairs(pcqn3,col=factor(targets$Array)),pch=20)

# cell type composition using houseman
cellCounts <-estimateCellCounts(rgSet,compositeCellType="Blood",cellTypes=c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),returnAll=T,meanPlot=TRUE,verbose=TRUE)
phenotype_data <- phenotype_data[order(match(rownames(phenotype_data),targets$Sample_ID)),]
cellCounts = cellCounts[[1]]
cellCounts <- cellCounts[order(match(rownames(cellCounts),targets$Sentrix_slide)),]

# batch effect correction on houseman data - not clear if necessary
cellCounts1 = t(ComBat(t(cellCounts),targets$Date_of_typing,mod0))
cellCounts2 = t(ComBat(t(cellCounts1),targets$Slide,mod0))
cellCounts3 = t(ComBat(t(cellCounts1),targets$Array,mod0))

plot(rowSums(ccomp),pch=20)
plot(rowSums(ccomp1),pch=20)

# scales sums of fractions to 100%
ccomp = t(apply(ccomp,1,function(x) x/sum(x)))

pairs(cellCounts1,col=factor(targets$Date_of_typing),pch=20)
pairs(cellCounts2,col=factor(targets$Slide),pch=20)
pairs(cellCounts3,col=factor(targets$Array),pch=20)

pairs(apply(cellCounts3,2,scale),col=factor(pheno0$Array),pch=20)

save(mfn,mqn,mfn1,mqn1,mfn2,mqn2,mfn3,mqn3,cellCounts,cellCounts1,cellCounts2,cellCounts3,file="media/sugi/NGS_Core_Bonn_1/Extremegroups/data/batch_corrected_cellcounts_extremegroups.Rdata",compress=FALSE)



