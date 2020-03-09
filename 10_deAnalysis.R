# File: 10_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 13/12/2019

## load the data
source('header.R')
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# select the right table using data and project id
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 43) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n[2])

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 43')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## make count matrix
names(lCounts)
mCounts = do.call(cbind, lCounts)
colnames(mCounts) = names(lCounts)

# sanity check
identical(dfSample$id, as.integer(colnames(mCounts)))

mData = mCounts
dim(mData)

## use combination of batch and biological source as identifier for technical replicates
fReplicates = factor(dfSample$group1):factor(dfSample$group2)
levels(fReplicates)
dfSample$fReplicates = factor(fReplicates)
# combine the technical replicates
i = seq_along(1:ncol(mData))
m = tapply(i, dfSample$fReplicates, function(x) {
  return(x)
})

mData = sapply(m, function(x){
  return(rowSums(mCounts[,x]))
})

# get a shorter version of dfSample after adding technical replicates
dfSample.2 = dfSample[sapply(m, function(x) return(x[1])), ]
identical(colnames(mData), as.character(dfSample.2$fReplicates))
dim(dfSample.2)
dfSample.2 = droplevels.data.frame(dfSample.2)

## normalise the data
# drop the rows where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
# FALSE  TRUE 
# 20254  6941 
mData = mData[!(i< 3),]
dim(mData)
# [1] 20254     4

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

identical(colnames(mData.norm), as.character(dfSample.2$fReplicates))

## delete sample section after testing
mData.norm = round(mData.norm, 0)

# set.seed(123)
# i = sample(1:nrow(mData.norm), 100, replace = F)
# dfData = data.frame(t(mData.norm[i,]))

dfData = data.frame(t(mData.norm))
dim(dfData)
dfData = stack(dfData)

## create covariates for modelling
str(dfSample.2)
dfData$fTreatment = factor(dfSample.2$group1)
dfData$fPatient = factor(dfSample.2$group2)
dfData = droplevels.data.frame(dfData)

dfData$Coef.1 = factor(dfData$fTreatment:dfData$ind)
dfData$Coef.2 = factor(dfData$fPatient:dfData$ind)
str(dfData)

# # setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1 + (1 | Coef.1) + (1 | Coef.2), data=dfData)
# summary(fit.lme1)
# fit.lme2 = glmer.nb(values ~ 1 + (1 | Coef.1), data=dfData)
# summary(fit.lme2)
# anova(fit.lme1, fit.lme2)
# ran = ranef(fit.lme1, condVar=F)
# 
# plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp2RandomEffectsMultipleScales.stan')

## calculate hyperparameters for variance of coefficients
# l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))
# # ## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# r2 = ran$Coef.adj1
# r3 = ran$Coef.adj2
# 
# initf = function(chain_id = 1) {
#   list(sigmaRan1 = 1, sigmaRan2=1)
# }

## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef.1), ]
d2 = dfData[!duplicated(dfData$Coef.2), ]

lStanData = list(Ntotal=nrow(dfData), 
                 Nclusters1=nlevels(dfData$Coef.1),
                 Nclusters2=nlevels(dfData$Coef.2),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NScaleBatches2 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef.1),
                 NgroupMap2=as.numeric(dfData$Coef.2),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 NBatchMap2=as.numeric(d2$ind), # this is where we use the second level mapping
                 Nphi=nlevels(dfData$ind),
                 NphiMap=as.numeric(dfData$ind),
                 y=dfData$values, 
                 #gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*3)

ptm = proc.time()

fit.stan = sampling(stanDso, data=lStanData, iter=1500, chains=4,
                    pars=c('sigmaRan1',
                           'sigmaRan2',
                           'phi',
                           #'mu',
                           'rGroupsJitter1',
                           'rGroupsJitter2',
                           'betas'
                           #'phi_scaled'
                           ),
                    cores=4, control=list(adapt_delta=0.99, max_treedepth = 11))#, init=initf)
save(fit.stan, file='results/fit.stan.nb_3Mar.rds')
ptm.end = proc.time()
print(fit.stan, c('sigmaRan1'), digits=3)
print(fit.stan, c('phi'), digits=3)
print(fit.stan, c('rGroupsJitter1'))
traceplot(fit.stan, c('sigmaRan1[1]'))
traceplot(fit.stan, c('sigmaRan1[2]'))
traceplot(fit.stan, c('rGroupsJitter1[1]', 'sigmaRan1[1]'))

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# # ## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan)$betas)
# ## add the intercept to each random effect variable, to get the full coefficient
# mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef.1))
# the split is done below on : symbol, but factor name has a : symbol due
# to creation of interaction earlier, do some acrobatics to sort that issue
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
str(d)
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='Non-lesional', deflection='Lesional') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
library(org.Hs.eg.db)
## remove X from annotation names
dfResults$ind = gsub('X', '', as.character(dfResults$ind))
df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'ENTREZID')
i = match(dfResults$ind, df$ENTREZID)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(dfResults$ind, df$ENTREZID)
## produce the plots 
f_plotVolcano(dfResults, 'lei vs nl')#, fc.lim=c(-2.5, 2.5))
f_plotVolcano(dfResults, 'lei vs nl', fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.01, 'lei vs nl')
table(dfResults$adj.P.Val < 0.01)
## save the results 
write.csv(dfResults, file='results/reverse/DEAnalysisLesionalVsNonLesional.xls')

######### do a comparison with deseq2
str(dfSample.2)
dfDesign = data.frame(Treatment = factor(dfSample.2$group1, levels = c('Non-lesional', 'Lesional')), Patient=factor(dfSample.2$group2),
                      row.names=colnames(mData))

oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment + Patient)
oDseq = DESeq(oDseq)

plotDispEsts(oDseq)
oRes = results(oDseq, contrast = c('Treatment', 'Lesional', 'Non-lesional'))
plotMA(oRes)
temp = as.data.frame(oRes)
i = match((dfResults$ind), rownames(temp))
temp = temp[i,]
identical((dfResults$ind), rownames(temp))
plot(dfResults$logFC, log(2^temp$log2FoldChange), pch=20)
table(oRes$padj < 0.01)
write.csv(oRes, file='results/DESeq.xls')

r1 = dfResults
r2 = oRes

r1 = r1[order(abs(r1$logFC), decreasing = T),]
head(r1)
r2$logFC = log(2^r2$log2FoldChange)
r2 = r2[order(abs(r2$logFC), decreasing = T),]

r1.top = r1$ind[1:100]
r2.top = rownames(r2)[1:100]
table(r1.top %in% r2.top)
