# linear regression analysis using limma
# with and without combat

# variables of interest
subgroup <- as.factor(targets$subgroupSubgruppe)
age <- as.numeric(targets$age)
currSmoker <- as.factor(targets$currSmoker)
earlySmoker <- as.factor(targets$earlySmoker)
antibabyPill <- as.factor(targets$antibabyPill)
orgDisease <- as.factor(targets$orgDisease)
medication <- as.factor(targets$otherMedication)
CD8T <-  as.numeric(cellCounts3$CD8T)
CD4T <- as.numeric(cellCounts3$CD4T)
NK <- as.numeric(cellCounts3$NK)
Bcell <- as.numeric(cellCounts$Bcell)
Mono <- as.numeric(cellCounts3$Mono)
ccomp <- cbind(CD8T,CD4T,NK,Bcell,Mono)

# design matrix (no intercept included)
design <- model.matrix(~0+subgroup,data=targets)
design2 <- model.matrix(~0+subgroup+age+currSmoker+earlySmoker+antibabyPill+orgDisease+medication, data=targets)
design3 <- model.matrix(~0+subgroup+age+currSmoker+earlySmoker+antibabyPill+orgDisease+medication+cccomp,data=targets)

# fit the linear model with combat for quantile and funnorm normalized data
fit <- lmFit(mfn3, design)
fit <- lmFit(mqn3, design)

# create contrast matrix for group comparisons
contMatrix <- makeContrasts(subgroupControl-subgroupEnvironment_risk,subgroupControl-subgroupGenetic_risk,subgroupEnvironment_risk-subgroupGenetic_risk, levels=design)

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
topTable(fit2,coef=1,adjust="BH")

# fit the linear model without combat for quantile and funnorm normalized data
fit_without_combat <- lmFit(mfn, design)
fit_without_combat <- lmFit(mqn, design)
fit_without_combat2 <- contrasts.fit(fit, contMatrix)
fit_without_combat2 <- eBayes(fit2)
topTable(fit_without_combat2,coef=1,adjust="BH")

# manual check if combat works properly (only funnorm)
# subsets adj. pvalues <0.1 from toptable
topfn <- subset(topTable(fit2),adj.P.Val<=0.10)
topfn_withou_combat <- subset(topTable(fit_without_combat2),adj.P.Val<=0.10)

# gets the methylation values of cpgs with a adj. pvalue <0.1
topmfn <- data.frame(t(mfn[row.names(mfn) %in% row.names(topfn),]))
topmfn_without_combat <- data.frame(t(mfn[row.names(mfn) %in% row.names(topfn),]))

# manual calculation of residuals over combat corrected and not corrected data 
mfnlm <- list()
mfnlm_without_combat <- list()

for (i in row.names(topfn)) mfnlm[[i]] = {
  cpg <- topmfn[[i]]
  residuals <- residuals(lm(cpg ~ subgroup+targets$Date_of_typing+targets$Slide+targets$Array,data=targets))
  summary(residuals)
}
for (i in row.names(topfn_withou_combat)) mfnlm_without_combat[[i]] = {
  cpg <- topmfn_without_combat[[i]]
  residuals_without_combat = residuals(lm(cpg ~ subgroup+targets$Date_of_typing+targets$Slide+targets$Array,data=targets))
  summary(residuals)
}



              

